"""Inputs for the source-attributed exposure model (3 November 2026).

  data/sources/fire_grid_{YYYY}_{MM}.tif   daily MODIS FIRMS fire-pixel counts on a 0.05 degree grid, one band per day,
                                           27.0 to 32.5 E, 21.9 to 16.4 S (province plus about 100 km), 2017 to 2022
  data/sources/acag_monthly_districts.csv  ACAG V5.GL monthly PM2.5, area mean per district, 2017 to 2022 (second product
                                           for the monthly fire model)
Usage: python analysis/33_model_inputs.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse, json, time
from concurrent.futures import ThreadPoolExecutor
import ee
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'sources'
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
ee.Initialize(project=args.project)
W, N, STEP, NX, NY = 27.0, -16.4, 0.05, 110, 110
FIRMS = ee.ImageCollection('FIRMS').select('T21')
grid_proj = ee.Projection('EPSG:4326').atScale(1000)

def retry(fn, what):
    for attempt in range(6):
        try: return fn()
        except Exception as err:
            print(what, 'attempt', attempt + 1, str(err)[:90], flush=True); time.sleep(20)
    raise SystemExit(f'{what}: Earth Engine unavailable, rerun to resume')

def month_grid(ym):
    y, m = ym; f = OUT / f'fire_grid_{y}_{m:02d}.tif'
    if f.exists(): return
    start = ee.Date.fromYMD(y, m, 1); n = start.advance(1, 'month').difference(start, 'day').getInfo()
    def day(k):
        d0 = start.advance(k, 'day')
        img = FIRMS.filterDate(d0, d0.advance(1, 'day')).max().gt(0).unmask(0)
        img = ee.Image(ee.Algorithms.If(FIRMS.filterDate(d0, d0.advance(1, 'day')).size().gt(0), img, ee.Image.constant(0)))
        return img.rename('f').toFloat().setDefaultProjection(grid_proj) \
                  .reduceResolution(ee.Reducer.sum(), maxPixels=65535).rename(d0.format('YYYY-MM-dd'))
    stack = ee.ImageCollection([day(k) for k in range(n)]).toBands()
    data = retry(lambda: ee.data.computePixels({
        'expression': stack.unmask(0), 'fileFormat': 'GEO_TIFF',
        'grid': {'dimensions': {'width': NX, 'height': NY}, 'crsCode': 'EPSG:4326',
                 'affineTransform': {'scaleX': STEP, 'shearX': 0, 'translateX': W, 'shearY': 0, 'scaleY': -STEP, 'translateY': N}}}),
        f'fires {y}-{m:02d}')
    f.write_bytes(data); print('fires', y, m, flush=True)

with ThreadPoolExecutor(2) as pool:  # noncommercial concurrency limit
    list(pool.map(month_grid, [(y, m) for y in range(2017, 2023) for m in range(1, 13)]))

f = OUT / 'acag_monthly_districts.csv'
if not f.exists():
    fc = ee.FeatureCollection(json.loads((ROOT / 'data/public/study_administrative_units.geojson').read_text()))
    col = ee.ImageCollection('projects/sat-io/open-datasets/GLOBAL-SATELLITE-PM25/MONTHLY').filterDate('2017-01-01', '2023-01-01')
    rows = []
    for y in range(2017, 2023):
        sub = col.filterDate(f'{y}-01-01', f'{y + 1}-01-01')
        def red(img):
            return img.reduceRegions(fc, ee.Reducer.mean(), 1000).map(
                lambda ft: ft.set('month', img.date().format('YYYY-MM')).select(['shapeName', 'month', 'mean']))
        res = retry(lambda: sub.map(red).flatten().getInfo()['features'], f'acag {y}')
        rows += [x['properties'] for x in res]; print('acag', y, flush=True)
    pd.DataFrame(rows).rename(columns={'shapeName': 'area', 'mean': 'pm25'}).to_csv(f, index=False)
print('inputs complete')
