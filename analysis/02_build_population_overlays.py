"""Prepare reusable native-grid population weights for all study years."""
from pathlib import Path
import json,os,sys
os.environ.pop('PROJ_LIB',None);os.environ.pop('PROJ_DATA',None)
os.environ['OPENBLAS_NUM_THREADS']='1'
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np,pandas as pd,rasterio
from shapely.geometry import shape
from shapely.ops import unary_union
from population_overlay import geometry,union_population_support,matrix_for_area,apply_counts
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data/public';QA=ROOT/'qa';OUT=ROOT/'results'
CACHE=QA/'overlay_cache';CACHE.mkdir(exist_ok=True)
cfg=json.loads((ROOT/'analysis/study_config.json').read_text())
features=json.loads((P/'study_administrative_units.geojson').read_text())['features']
boundaries={f['properties']['shapeName']:shape(f['geometry']) for f in features}
names=list(boundaries);TOTAL='Midlands Province';boundaries[TOTAL]=unary_union(list(boundaries.values()))
def key(name):return name.replace(' ','_')
paths=[P/f'worldpop_population_{y}.tif' for y in cfg['years']]
assert all(p.exists() for p in paths),'All year-specific population inputs must be downloaded first'
support=union_population_support(paths,CACHE)
print('Positive fine-grid support across all years:',len(support['rows']),flush=True)
summaries=[];totals={}
for product,prefix in [('primary','ghap'),('comparator','acag')]:
    with rasterio.open(P/f'{prefix}_annual_2022.tif') as src:profile=src.profile.copy()
    geos=geometry(profile,boundaries)
    area_sum=sum(geos[n]['area_m2'] for n in names)
    assert np.max(np.abs(area_sum-geos[TOTAL]['area_m2']))<.02,'Reporting areas must partition provincial cell footprints'
    np.savez_compressed(QA/f'spatial_support_{product}.npz',area_names=np.array(list(boundaries)),
                        **{key(n):g['area_m2'] for n,g in geos.items()})
    matrices={}
    for name in names:
        matrices[name]=matrix_for_area(support,geos[name],profile,boundaries[name],CACHE,product+'_'+name)
    for year in cfg['years']:
        with rasterio.open(P/f'worldpop_population_{year}.tif') as src:pop=src.read(1).astype('float64')
        weights={};centre={}
        for name in names:
            weights[name],centre[name],qa=apply_counts(matrices[name],pop)
            summaries.append({'year':year,'area':name,'product':product,**qa})
        weights[TOTAL]=sum(weights.values());centre[TOTAL]=sum(centre.values())
        np.savez_compressed(OUT/f'population_weights_{product}_{year}.npz',
                            **{key(n):v for n,v in weights.items()},
                            **{'centre_'+key(n):v for n,v in centre.items()})
        total=float(weights[TOTAL].sum())
        if product=='primary':totals[year]=total
        else:assert abs(total-totals[year])<.01,(year,total,totals[year])
        print(product,year,'province population',round(total),flush=True)
    del geos,matrices
pd.DataFrame(summaries).to_csv(QA/'population_conservation.csv',index=False)
(QA/'population_overlay_summary.json').write_text(json.dumps({'positive_source_cells_all_years':len(support['rows']),
    'population_source_shape':[support['profile']['height'],support['profile']['width']],
    'method':'Sparse native-footprint intersection fractions in EPSG6933, applied to year-specific population counts',
    'population_by_year':totals,'district_names':names,'province_aggregation':'Sum of eight non-overlapping reporting areas'},indent=2))
print('All population overlays complete',flush=True)
