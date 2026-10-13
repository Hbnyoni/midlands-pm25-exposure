"""Download native-grid provincial subsets, with bounded tiles and full provenance.

Uses authenticated read-only Earth Engine access. No Drive/cloud exports, health
data uploads, interpolation, personal credentials or signed URLs are recorded.
"""
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,calendar,hashlib,json,math,os,threading,time,sys
os.environ.pop('PROJ_LIB',None);os.environ.pop('PROJ_DATA',None)
sys.stdout.reconfigure(encoding='utf-8')
import ee,requests,numpy as np,rasterio
from rasterio.io import MemoryFile
from rasterio.transform import Affine
from rasterio.windows import Window
from shapely.geometry import shape

ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data/public';M=ROOT/'data/metadata'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--project',required=True)
parser.add_argument('--year',type=int,action='append')
parser.add_argument('--month',type=int,action='append')
parser.add_argument('--workers',type=int,default=3)
args=parser.parse_args()
cfg=json.loads((ROOT/'analysis/study_config.json').read_text())
ee.Initialize(project=args.project)
province=shape(json.loads((P/'midlands_province.geojson').read_text())['features'][0]['geometry'])
x0,y0,x1,y1=province.bounds
extent=[x0-.025,y0-.025,x1+.025,y1+.025]
collections={'daily':ee.ImageCollection(cfg['ghap_daily_collection']),
             'annual':ee.ImageCollection(cfg['ghap_annual_collection']),
             'population':ee.ImageCollection(cfg['worldpop_collection']).filter(ee.Filter.eq('country','ZWE')).filter(ee.Filter.eq('release','R2025A')),
             'comparator':ee.ImageCollection('projects/sat-io/open-datasets/GLOBAL-SATELLITE-PM25/ANNUAL')}
catalog={}
for label,col in collections.items():
    first=col.first();catalog[label]={'first_image':first.getInfo(),'projection':first.select([0]).projection().getInfo()}
    print(label,catalog[label]['first_image']['id'],flush=True)
(M/'gee_collection_discovery.json').write_text(json.dumps({'sources':catalog,'extent':extent,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'survey_inputs':False},indent=2),encoding='utf-8')
manifest_path=M/'gee_download_manifest.json'
manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else []
lock=threading.Lock()
api_lock=threading.Lock()  # Earth Engine's shared discovery transport is serial.
def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def fetch(job):
    image,destination,ids,projection,details=job
    rel=str(destination.relative_to(ROOT))
    old=next((z for z in manifest if z['file']==rel),None)
    if old and all(old.get(k)==v for k,v in details.items()) and old['source_asset_ids']==ids and old['region_bounds_wgs84']==extent and destination.exists() and digest(destination)==old['sha256']:
        print('Verified',destination.name,flush=True);return
    t=Affine(*projection['transform']);assert t.b==0 and t.d==0 and t.a>0 and t.e<0
    c0=math.floor((extent[0]-t.c)/t.a);c1=math.ceil((extent[2]-t.c)/t.a)
    r0=math.floor((extent[3]-t.f)/t.e);r1=math.ceil((extent[1]-t.f)/t.e)
    width,height=c1-c0,r1-r0;transform=t*Affine.translation(c0,r0)
    bands=len(details.get('band_dates',[None]))
    edge=min(2000,int(math.sqrt(22_000_000/(4*bands))))
    temp=destination.with_suffix('.partial.tif');tiles=[]
    profile={'driver':'GTiff','height':height,'width':width,'count':bands,'dtype':'float32','crs':projection['crs'],
             'transform':transform,'nodata':-9999,'compress':'deflate','predictor':3,'tiled':True}
    with rasterio.open(temp,'w',**profile) as output:
        for row in range(0,height,edge):
            for col in range(0,width,edge):
                w,h=min(edge,width-col),min(edge,height-row)
                tile_transform=transform*Affine.translation(col,row)
                params={'crs':projection['crs'],'crs_transform':list(tile_transform)[:6],
                        'dimensions':[w,h],'format':'GEO_TIFF','filePerBand':False}
                for attempt in range(6):
                    try:
                        with api_lock:
                            url=image.toFloat().unmask(-9999).getDownloadURL(params)
                        response=requests.get(url,timeout=(30,240));response.raise_for_status()
                        with MemoryFile(response.content) as memory:
                            with memory.open() as tile:
                                assert (tile.width,tile.height,tile.count)==(w,h,bands),(tile.shape,w,h,bands)
                                assert tile.transform.almost_equals(tile_transform,precision=1e-9),(tile.transform,tile_transform)
                                output.write(tile.read(),window=Window(col,row,w,h))
                        tiles.append({'row':row,'column':col,'width':w,'height':h,'download_sha256':hashlib.sha256(response.content).hexdigest()})
                        break
                    except (requests.RequestException,ee.EEException,ConnectionError,TimeoutError) as error:
                        if attempt==5:raise
                        print('Retry',destination.name,type(error).__name__,flush=True)
                        time.sleep(2**attempt)
        output.update_tags(source_asset_ids=json.dumps(ids),product=details['product'],units=details['units'],native_grid='preserved without interpolation')
        if 'band_dates' in details:
            for b,date in enumerate(details['band_dates'],1):output.set_band_description(b,date)
    temp.replace(destination)
    record={'file':rel,'sha256':digest(destination),'bytes':destination.stat().st_size,'source_asset_ids':ids,
            'projection':projection,'extract_transform':list(transform)[:6],'extract_dimensions':[width,height],
            'region_bounds_wgs84':extent,'retrieved_utc':datetime.now(timezone.utc).isoformat(),
            'tile_assembly':'Integer windows on unchanged native grid; lossless compression; no resampling','tiles':tiles,**details}
    with lock:
        manifest[:]=[m for m in manifest if m['file']!=rel]+[record]
        scratch=manifest_path.with_suffix('.new.json');scratch.write_text(json.dumps(manifest,indent=2),encoding='utf-8');scratch.replace(manifest_path)
    print('Downloaded',destination.name,width,height,f'{destination.stat().st_size/1e6:.1f} MB',flush=True)

jobs=[]
for year in (args.year or cfg['years']):
    start=f'{year}-01-01';end=f'{year+1}-01-01'
    sub=collections['daily'].filterDate(start,end).sort('system:time_start')
    ids=sub.aggregate_array('system:index').getInfo();times=sub.aggregate_array('system:time_start').getInfo()
    dates=[datetime.fromtimestamp(t/1000,timezone.utc).strftime('%Y-%m-%d') for t in times]
    assert len(dates)==len(set(dates))
    (M/f'ghap_daily_index_{year}.json').write_text(json.dumps({'year':year,'expected_days':366 if calendar.isleap(year) else 365,'available_days':len(dates),'images':[{'date':d,'id':i} for d,i in zip(dates,ids)]},indent=2))
    for month in (args.month or range(1,13)):
        lo=f'{year}-{month:02d}-01';hi=f'{year+1}-01-01' if month==12 else f'{year}-{month+1:02d}-01'
        entries=[(d,i) for d,i in zip(dates,ids) if lo<=d<hi]
        image=sub.filterDate(lo,hi).sort('system:time_start').select([0]).toBands().multiply(cfg['ghap_scale_factor'])
        jobs.append((image,P/f'ghap_daily_{year}_{month:02d}.tif',[i for d,i in entries],catalog['daily']['projection'],
                     {'product':'GHAP daily surface PM2.5','band_dates':[d for d,i in entries],'units':'micrograms/m3','source_scale_factor_applied':cfg['ghap_scale_factor'],'nodata':-9999}))
    for key,prefix,product,units in [('annual','ghap_annual','GHAP annual PM2.5','micrograms/m3'),('population','worldpop_population','WorldPop public population count','persons per native cell'),('comparator','acag_annual','ACAG annual surface PM2.5 comparison','micrograms/m3')]:
        col=collections[key].filterDate(start,end);items=col.aggregate_array('system:index').getInfo()
        assert len(items)==1,(key,year,items)
        image=col.first().select('population' if key=='population' else [0]);pr=image.projection().getInfo()
        detail={'product':product,'year':year,'units':units,'nodata':-9999}
        if key=='annual':detail['source_scale_factor_applied']=cfg['ghap_scale_factor'];image=image.multiply(cfg['ghap_scale_factor'])
        if key=='population':
            assert 'ZWE_' in items[0] and '_CN_' in items[0] and '_R2025A_V1' in items[0],items
            detail['release']='R2025A v1 constrained; Zimbabwe'
        jobs.append((image,P/f'{prefix}_{year}.tif',items,pr,detail))

jobs.sort(key=lambda j:0 if j[4]['product'].startswith('WorldPop') else 1 if 'annual' in j[4]['product'] else 2)
print('Prepared',len(jobs),'independent raster downloads',flush=True)
with ThreadPoolExecutor(max_workers=max(1,min(args.workers,4))) as pool:
    futures=[pool.submit(fetch,job) for job in jobs]
    for future in as_completed(futures):future.result()
print('Province acquisition complete',flush=True)
