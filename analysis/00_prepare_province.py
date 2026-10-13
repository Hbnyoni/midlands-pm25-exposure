"""Prepare an independently sourced Midlands Province study geography.

The Shurugwi public-data study is retained as a separate preliminary analysis.
No survey or simulated dataset is accessed.
"""
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib, shutil, os, re
os.environ.pop('PROJ_LIB',None);os.environ.pop('PROJ_DATA',None)
import requests, shapely
from shapely.geometry import shape, mapping
from shapely.ops import unary_union, transform
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[1]
for part in ['data/public','data/metadata','results/rasters','qa','figures','evidence','manuscript']:
    (ROOT/part).mkdir(parents=True,exist_ok=True)
P=ROOT/'data/public';M=ROOT/'data/metadata'
# The archive supplies its own helpers, configuration and evidence records.
# Reproduction does not depend on the retained preliminary Shurugwi package.
for relative in ['evidence/published_risk_functions.csv','data/metadata/ghap_units_audit.json','data/public/natural_earth_110m_countries.geojson','analysis/study_config.json']:
    if not (ROOT/relative).exists():raise FileNotFoundError(f'Missing packaged prerequisite: {relative}')
config=json.loads((ROOT/'analysis/study_config.json').read_text())
config.update({'design':'Historical population exposure assessment with supplementary literature-transferred mortality scenarios',
    'working_geography':'Midlands Province, Zimbabwe; public administrative subunits compared within the province',
    'geography_status':'Province expansion adopted following user instruction to continue the recommended redesign on 2026-09-23',
    'primary_results':['Annual and daily WHO exposure comparisons','Population and population-day distributions','Administrative subunit comparison','Exposure-product sensitivity'],
    'mortality_analysis_role':'Supplementary scenarios only; no observed mortality endpoint',
    'figure_text':'Lowercase panel letters only; descriptive titles and source notes in manuscript captions',
    'analysis_plan_date':'2026-09-23'})
sources=[]
def fetch(level):
    url=f'https://raw.githubusercontent.com/wmgeolab/geoBoundaries/9469f09/releaseData/gbOpen/ZWE/{level}/geoBoundaries-ZWE-{level}.geojson'
    path=P/f'zimbabwe_{level.lower()}.geojson'
    if not path.exists():
        response=requests.get(url,timeout=(30,180));response.raise_for_status()
        temp=path.with_suffix('.partial');temp.write_bytes(response.content);temp.replace(path)
    if path.read_bytes()[:50].startswith(b'version https://git-lfs.github.com/spec/v1'):
        pointer=path.read_text()
        expected=re.search(r'oid sha256:([0-9a-f]{64})',pointer).group(1)
        media=url.replace('raw.githubusercontent.com','media.githubusercontent.com/media')
        response=requests.get(media,timeout=(30,180));response.raise_for_status()
        assert hashlib.sha256(response.content).hexdigest()==expected,'Git LFS content hash mismatch'
        temp=path.with_suffix('.partial');temp.write_bytes(response.content);temp.replace(path)
    payload=json.loads(path.read_text(encoding='utf-8'))
    sources.append({'level':level,'source_url':url,'file':str(path.relative_to(ROOT)),
                    'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size,
                    'represented_year':2020,'source':'ZIMSTAT/OCHA via geoBoundaries gbOpen',
                    'licence':'CC BY 3.0 IGO','retrieved_utc':datetime.now(timezone.utc).isoformat()})
    return payload
provinces=fetch('ADM1');districts=fetch('ADM2')
chosen=[f for f in provinces['features'] if f['properties']['shapeName'].casefold()=='midlands']
assert len(chosen)==1
province=shape(chosen[0]['geometry']);assert province.is_valid
tx=Transformer.from_crs(4326,6933,always_xy=True).transform
province_equal=transform(tx,province)
selected=[];records=[]
for f in districts['features']:
    g=shape(f['geometry'])
    if not g.intersects(province):continue
    eg=transform(tx,g);overlap=eg.intersection(province_equal)
    frac=overlap.area/eg.area
    records.append({'name':f['properties']['shapeName'],'fraction_of_adm2_in_province':frac,'overlap_km2':overlap.area/1e6})
    if frac>0.5:
        clipped=g.intersection(province)
        selected.append({'type':'Feature','properties':{**f['properties'],'province':'Midlands','province_overlap_fraction':frac},'geometry':mapping(clipped)})
assert selected
union=unary_union([shape(f['geometry']) for f in selected])
gap=province.difference(union);overlaps=sum(transform(tx,shape(f['geometry'])).area for f in selected)-transform(tx,union).area
gap_area=transform(tx,gap).area
# Source hierarchy must partition the target. Do not silently discard gaps.
if gap_area>province_equal.area*1e-5 or overlaps>province_equal.area*1e-5:
    raise RuntimeError(f'Administrative topology needs review: gap {gap_area/1e6} km2; overlaps {overlaps/1e6} km2')
groups={
    'Chirumhanzu':['Chirumhanzu'],
    'Gokwe North':['Gokwe North'],
    'Gokwe South':['Gokwe South','Gokwe South Urban'],
    'Gweru':['Gweru','Gweru Urban'],
    'Kwekwe':['Kwekwe','Kwekwe Urban','Redcliff'],
    'Mberengwa':['Mberengwa'],
    'Shurugwi':['Shurugwi','Shurugwi Town'],
    'Zvishavane':['Zvishavane','Zvishavane Urban']}
assert sorted(n for names in groups.values() for n in names)==sorted(f['properties']['shapeName'] for f in selected)
source_units={f['properties']['shapeName']:shape(f['geometry']) for f in selected}
grouped=[{'type':'Feature','properties':{'shapeName':name,'province':'Midlands','contributing_source_units':names,
           'geography_role':'District reporting area with associated municipal units combined'},
          'geometry':mapping(unary_union([source_units[n] for n in names]))} for name,names in groups.items()]
study={'type':'FeatureCollection','features':grouped}
(P/'source_adm2_units.geojson').write_text(json.dumps({'type':'FeatureCollection','features':selected}),encoding='utf-8')
(P/'study_administrative_units.geojson').write_text(json.dumps(study),encoding='utf-8')
(P/'midlands_province.geojson').write_text(json.dumps({'type':'FeatureCollection','features':chosen}),encoding='utf-8')
config['boundary_names']=list(groups)
config['district_reporting_groups']=groups
(ROOT/'analysis/study_config.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
record={'sources':sources,'province':'Midlands','district_selection':'More than half of source ADM2 equal-area polygon lies inside Midlands; selected polygons clipped to province',
    'selected_names':config['boundary_names'],'source_adm2_names':list(source_units),'district_reporting_groups':groups,
    'grouping_note':'Eight district reporting areas combine rural and associated municipal source units; these analytical merges are explicitly distinguished from the 14 unmodified source ADM2 units',
    'all_intersections':records,'province_area_km2':province_equal.area/1e6,
    'province_bounds':list(province.bounds),'unassigned_gap_km2':gap_area/1e6,'overlapping_adm2_area_km2':overlaps/1e6,
    'survey_coordinates_used':False,'subset_sha256':hashlib.sha256((P/'study_administrative_units.geojson').read_bytes()).hexdigest()}
(M/'boundary_provenance.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps({k:record[k] for k in ['province_area_km2','province_bounds','selected_names','unassigned_gap_km2','overlapping_adm2_area_km2']},indent=2),flush=True)
