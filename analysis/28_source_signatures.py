"""Source signatures behind modelled PM2.5 exposure in Midlands Province (3 November 2026, GeoHealth resubmission).

Earth Engine acquisition (cached in data/sources/):
  tropomi_{no2,so2}_{pre,post}.tif   TROPOMI OFFL L3 tropospheric NO2 and total SO2 columns, mean July 2018 to May 2024 (pre)
                                     and July 2024 to September 2026 (post), 0.03 degree grid over Midlands + 1.5 degree margin
  fires_daily.csv                    MODIS FIRMS fire pixels per day within Midlands + 100 km, 2017 to 2022
  era5_daily.csv                     ERA5-Land daily total precipitation and 10 m wind speed, Midlands mean, 2017 to 2022
  monthly_sites.csv                  monthly NO2 and SO2 means within 10 km of each detected hotspot and in a 20 to 50 km ring,
                                     July 2018 to September 2026
Sites are taken at their OpenStreetMap coordinates (Manhize, Redcliff ZISCO, Unki, Mimosa) plus Gweru town; Kwekwe is the
largest pre-period NO2 column within 25 km of the town (Zimasco is not mapped). Usage: python analysis/28_source_signatures.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse, json, time
import ee
import numpy as np
import pandas as pd
import rasterio

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'sources'; OUT.mkdir(parents=True, exist_ok=True)
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
ee.Initialize(project=args.project)

prov = json.load(open(ROOT / 'data/public/midlands_province.geojson'))
PROV = ee.Geometry(prov['features'][0]['geometry'])
REG = PROV.buffer(100000)
W, S, E, N = 27.5, -21.8, 32.0, -17.2          # Midlands plus a margin that holds every 50 km ring
STEP = 0.03
NO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').select('tropospheric_NO2_column_number_density')
SO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_SO2').select('SO2_column_number_density')
PRE, POST = ('2018-07-01', '2024-06-01'), ('2024-07-01', '2026-10-01')   # Manhize blast furnace fired 12 June 2024

def retry(fn, what):
    for attempt in range(5):
        try: return fn()
        except Exception as err:
            print(what, 'attempt', attempt + 1, str(err)[:90], flush=True); time.sleep(20)
    raise SystemExit(f'{what}: Earth Engine unavailable, rerun to resume')

def grid(name, img):
    f = OUT / f'{name}.tif'
    if f.exists(): return f
    nx, ny = int(round((E - W) / STEP)), int(round((N - S) / STEP))
    data = retry(lambda: ee.data.computePixels({
        'expression': img.multiply(1e6).rename('v').toFloat().unmask(-9999), 'fileFormat': 'GEO_TIFF',
        'grid': {'dimensions': {'width': nx, 'height': ny}, 'crsCode': 'EPSG:4326',
                 'affineTransform': {'scaleX': STEP, 'shearX': 0, 'translateX': W, 'shearY': 0, 'scaleY': -STEP, 'translateY': N}}}), name)
    f.write_bytes(data); print(name, flush=True); return f

for gas, col in (('no2', NO2), ('so2', SO2)):
    grid(f'tropomi_{gas}_pre', col.filterDate(*PRE).mean())
    grid(f'tropomi_{gas}_post', col.filterDate(*POST).mean())

def rd(f):
    with rasterio.open(f) as src:
        a = src.read(1).astype(float); t = src.transform
    a[a < -9000] = np.nan
    return a, t

# Sites at mapped coordinates (OpenStreetMap, checked 3 November 2026). Zimasco Kwekwe is not mapped, so the Kwekwe point is
# the largest pre-period NO2 column within 25 km of the town centre. A data-driven 'largest increase' search was rejected
# because it selected Unki rather than Manhize, which lies just outside the province in Chikomba district.
pre, t = rd(OUT / 'tropomi_no2_pre.tif')
ys, xs = np.mgrid[0:pre.shape[0], 0:pre.shape[1]]
lon = t.c + (xs + 0.5) * t.a; lat = t.f + (ys + 0.5) * t.e
kw = np.where(np.hypot((lon - 29.80) * 111.32 * np.cos(np.radians(-18.95)), (lat + 18.95) * 110.57) <= 25, pre, np.nan)
iy, ix = np.unravel_index(np.nanargmax(kw), pre.shape)
sites = {'Manhize': (30.6085, -18.9441), 'Kwekwe': (round(float(lon[iy, ix]), 4), round(float(lat[iy, ix]), 4)),
         'Redcliff': (29.7503, -19.0157), 'Unki': (30.0955, -19.6207), 'Mimosa': (29.8215, -20.3271), 'Gweru': (29.82, -19.45)}
json.dump(sites, open(OUT / 'hotspots.json', 'w'), indent=1); print('hotspots', sites)

f = OUT / 'monthly_sites.csv'
if not f.exists():
    others = ee.Geometry.MultiPoint([list(v) for v in sites.values()]).buffer(15000)
    rows = []
    months = pd.date_range('2018-07-01', '2026-09-01', freq='MS')
    for name, (lo, la) in sites.items():
        pt = ee.Geometry.Point(lo, la)
        near, ring = pt.buffer(10000), pt.buffer(50000).difference(pt.buffer(20000)).difference(others)
        def month_stats(m):
            m = ee.Date(m)
            empty = lambda b: ee.ImageCollection([ee.Image.constant(0).rename(b).toFloat().updateMask(0)])
            # months with no valid retrievals return a fully masked band, so the reduction gives null instead of failing
            img = ee.Image.cat([NO2.filterDate(m, m.advance(1, 'month')).merge(empty('tropospheric_NO2_column_number_density')).mean().rename('no2'),
                                SO2.filterDate(m, m.advance(1, 'month')).merge(empty('SO2_column_number_density')).mean().rename('so2')]).multiply(1e6)
            a = img.reduceRegion(ee.Reducer.mean(), near, 3500); b = img.reduceRegion(ee.Reducer.mean(), ring, 3500)
            return ee.Feature(None, {'month': m.format('YYYY-MM'), 'no2_near': a.get('no2'), 'so2_near': a.get('so2'),
                                     'no2_ring': b.get('no2'), 'so2_ring': b.get('so2')})
        fc = ee.FeatureCollection(ee.List([str(m.date()) for m in months]).map(month_stats))
        res = retry(lambda: fc.getInfo()['features'], f'monthly {name}')
        rows += [dict(x['properties'], site=name) for x in res]
        print('monthly', name, flush=True)
    pd.DataFrame(rows).to_csv(f, index=False)

f = OUT / 'fires_daily.csv'
if not f.exists():
    rows = []
    for y in range(2017, 2023):
        col = ee.ImageCollection('FIRMS').filterDate(f'{y}-01-01', f'{y + 1}-01-01').filterBounds(REG).select('T21')
        fc = col.map(lambda i: ee.Feature(None, {'date': i.date().format('YYYY-MM-dd'),
                                                 'fires': i.gt(0).reduceRegion(ee.Reducer.sum(), REG, 1000, maxPixels=1e9).get('T21')}))
        rows += [x['properties'] for x in retry(lambda: fc.getInfo()['features'], f'fires {y}')]
        print('fires', y, flush=True)
    pd.DataFrame(rows).groupby('date', as_index=False).fires.sum().to_csv(f, index=False)

f = OUT / 'era5_daily.csv'
if not f.exists():
    rows = []
    for y in range(2017, 2023):
        col = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').filterDate(f'{y}-01-01', f'{y + 1}-01-01')
        def stats(i):
            ws = i.select('u_component_of_wind_10m').hypot(i.select('v_component_of_wind_10m')).rename('wind')
            v = ee.Image.cat([i.select('total_precipitation_sum').multiply(1000).rename('precip_mm'), ws]) \
                .reduceRegion(ee.Reducer.mean(), PROV, 11000)
            return ee.Feature(None, {'date': i.date().format('YYYY-MM-dd'), 'precip_mm': v.get('precip_mm'), 'wind': v.get('wind')})
        rows += [x['properties'] for x in retry(lambda: col.map(stats).getInfo()['features'], f'era5 {y}')]
        print('era5', y, flush=True)
    pd.DataFrame(rows).to_csv(f, index=False)
print('acquisition complete')
