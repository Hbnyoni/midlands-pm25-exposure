"""Population-weighted spatial distributions, without smoothing or pseudo-replication.

An additional descriptive analysis after the initial manuscript. Quantiles are
the inverse weighted empirical CDF, not temporal quantiles or confidence limits.
"""
from pathlib import Path
from datetime import datetime, timezone
import os, json, shutil
os.environ.pop('PROJ_LIB', None); os.environ.pop('PROJ_DATA', None)
import numpy as np
import pandas as pd
import rasterio

ROOT=Path(__file__).resolve().parents[1]; RES=ROOT/'results'; QA=ROOT/'qa'
backup=ROOT/'revisions/before_graphics_revision_2026_09_30'
backup.mkdir(parents=True,exist_ok=True)
for name in ['Frontiers_Midlands_Review_Package.zip','analysis/figure_style.py',
             'analysis/07_build_manuscript.py','analysis/quality_manuscript.py']:
    src=ROOT/name; dst=backup/Path(name).name
    if src.exists() and not dst.exists(): shutil.copy2(src,dst)
plan=dict(recorded_utc=datetime.now(timezone.utc).isoformat(),
          status='Additional descriptive analysis after initial results; not preregistered',
          question='How much within-area spatial exposure variation is hidden by population-weighted means?',
          population='Year-matched, fractionally allocated WorldPop counts on each native exposure grid',
          rules=['Retain both annual products, all six years, province and all eight reporting areas',
                 'F(c) = sum(P_i * I(C_i <= c)) / sum(P_i)',
                 'Q(p) is the smallest native concentration for which F(c) >= p',
                 'Report spatial percentiles 10, 25, 50, 75 and 90 without kernel smoothing',
                 'Reconcile means, counts and WHO-threshold shares against the existing analysis',
                 'Spatial percentile ranges are descriptive exposure spread, never confidence intervals',
                 'No inferential tests on model cells and no new health coefficients'])
(ROOT/'evidence/distribution_analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
cells=pd.read_csv(RES/'native_cell_annual_metrics.csv',usecols=['year','cell_index','annual_pm25'])
baseline=pd.read_csv(RES/'annual_exposure_summary.csv')
comparison=pd.read_csv(RES/'exposure_product_sensitivity.csv')
names=baseline.area.drop_duplicates().tolist(); result=[]; curves=[]; raw=[]
for year in range(2017,2023):
    for product,kind in [('GHAP','primary'),('ACAG','comparator')]:
        with np.load(RES/f'population_weights_{kind}_{year}.npz') as z:
            weights={area:z[area.replace(' ','_')] for area in names}
        if product=='GHAP':
            values=np.full(len(weights[names[0]]),np.nan)
            g=cells[cells.year==year]
            values[g.cell_index.to_numpy()]=g.annual_pm25.to_numpy()
        else:
            with rasterio.open(ROOT/f'data/public/acag_annual_{year}.tif') as src:
                values=src.read(1,masked=True).filled(np.nan).ravel().astype(float)
        for area in names:
            w=weights[area]; valid=(w>0)&np.isfinite(values)
            assert np.isclose(w[valid].sum(),w.sum(),rtol=1e-12), (product,year,area)
            x=values[valid]; p=w[valid]; sort=np.argsort(x,kind='stable'); x=x[sort]; p=p[sort]
            total=p.sum(); c=np.cumsum(p)/total
            record=dict(year=year,product=product,area=area,population=total,
                        mean_pm25=np.sum(x*p)/total)
            for prob in [10,25,50,75,90]:
                record[f'p{prob}']=x[np.searchsorted(c,prob/100,side='left')]
            for cut in [5,10,15,25,35]:record[f'population_percent_above_{cut}']=100*p[x>cut].sum()/total
            record['p90_minus_p10']=record['p90']-record['p10']; result.append(record)
            old=baseline[(baseline.year==year)&(baseline.area==area)].iloc[0]
            oldcomp=comparison[(comparison.year==year)&(comparison.area==area)].iloc[0]
            expected=old.population_weighted_pm25 if product=='GHAP' else oldcomp.comparator_pw_pm25
            assert abs(record['mean_pm25']-expected)<1e-8
            assert abs(total-old.population)<1e-6
            expected_share=100*(old.population_above_annual_10 if product=='GHAP' else oldcomp.comparator_population_above_10)/old.population
            assert abs(record['population_percent_above_10']-expected_share)<1e-7
            # Plain native values and allocated counts permit independent base-R checks.
            raw.append(pd.DataFrame(dict(year=year,product=product,area=area,concentration=x,population=p)))
            if area=='Midlands Province':
                ends=np.r_[np.flatnonzero(np.diff(x)!=0),len(x)-1]
                curves.append(pd.DataFrame(dict(year=year,product=product,concentration=x[ends],population_at_or_below_percent=100*c[ends])))
pd.DataFrame(result).to_csv(RES/'population_exposure_distribution_summary.csv',index=False)
pd.concat(curves,ignore_index=True).to_csv(RES/'province_exposure_ecdf.csv',index=False)
pd.concat(raw,ignore_index=True).to_csv(QA/'distribution_native_inputs.csv',index=False)
latest=[r for r in result if r['year']==2022 and r['area']=='Midlands Province']
(RES/'distribution_results.json').write_text(json.dumps(dict(area_product_years=len(result),
    quantile_definition='Inverse population-weighted empirical CDF',province_2022=latest),indent=2),encoding='utf-8')
print(json.dumps(dict(area_product_years=len(result),province_2022=latest),indent=2))
