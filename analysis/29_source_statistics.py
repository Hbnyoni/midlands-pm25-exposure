"""Statistics for the source-signature analysis (3 November 2026, GeoHealth resubmission). Reads data/sources/ from 28.

1. Manhize steel plant (blast furnace fired 12 June 2024): monthly NO2 and SO2 excess (10 km disc minus 20 to 50 km ring),
   interrupted time series with calendar-month fixed effects and Newey-West errors (6 lags), June 2024 excluded as transition,
   at Manhize and at five Midlands sites (Kwekwe, Redcliff, Unki, Mimosa, Gweru). SO2 is reported but not interpreted: the
   regional SO2 background fell by about two thirds between periods, consistent with a retrieval change rather than emissions.
2. Fires and the September to November peak: daily provincial population-weighted GHAP PM2.5 (2017 to 2022) regressed on log
   fire pixels over the previous three days with year-month fixed effects, precipitation and wind (Newey-West, 7 lags), so the
   association is identified within months rather than from the shared seasonal cycle.
3. Kwekwe industrial belt: GHAP annual PM2.5 within 10 km of the NO2 hotspot against a 20 to 50 km ring, 2017 to 2022.
4. People near Manhize: WorldPop R2025A 2025 within 10 and 25 km of the hotspot.
Usage: python analysis/29_source_statistics.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
import rasterio
import statsmodels.formula.api as smf
import ee

ROOT = Path(__file__).resolve().parents[1]
SRC, RES = ROOT / 'data' / 'sources', ROOT / 'results'
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
sites = json.load(open(SRC / 'hotspots.json'))
OUT = {'hotspots': sites}

# ---- 1. Manhize interrupted time series -------------------------------------------------------------------------
m = pd.read_csv(SRC / 'monthly_sites.csv')
m['date'] = pd.to_datetime(m.month + '-01'); m['cal'] = m.date.dt.month
for g in ('no2', 'so2'): m[f'{g}_ex'] = m[f'{g}_near'] - m[f'{g}_ring']
m = m[m.date != '2024-06-01']
m['post'] = (m.date >= '2024-07-01').astype(int)
its = []
for site in sites:
    s = m[m.site == site].dropna(subset=['no2_ex'])
    for g in ('no2', 'so2'):
        d = s.dropna(subset=[f'{g}_ex'])
        fit = smf.ols(f'{g}_ex ~ post + C(cal)', d).fit(cov_type='HAC', cov_kwds={'maxlags': 6})
        pre_mean = d[d.post == 0][f'{g}_ex'].mean()
        lo, hi = fit.conf_int().loc['post']
        its.append({'site': site, 'gas': g, 'n_months': len(d), 'pre_mean': pre_mean, 'post_mean': d[d.post == 1][f'{g}_ex'].mean(),
                    'effect': fit.params['post'], 'ci_low': lo, 'ci_high': hi, 'p': fit.pvalues['post'],
                    'pct_of_pre': 100 * fit.params['post'] / abs(pre_mean) if pre_mean else np.nan})
its = pd.DataFrame(its); its.to_csv(RES / 'source_its.csv', index=False)
OUT['its'] = its.to_dict('records')
print(its.round(3).to_string(index=False))

# ---- 2. Fires and daily PM2.5 ----------------------------------------------------------------------------------
pm = pd.read_csv(RES / 'daily_area_exposure_and_scenarios.csv', parse_dates=['date'])
pm = pm[pm.area == 'Midlands Province'][['date', 'pw_pm25', 'population_days_above_15']]
fi = pd.read_csv(SRC / 'fires_daily.csv', parse_dates=['date'])
er = pd.read_csv(SRC / 'era5_daily.csv', parse_dates=['date'])
cal = pd.DataFrame({'date': pd.date_range('2017-01-01', '2022-12-31')})
d = cal.merge(fi, on='date', how='left').fillna({'fires': 0}).merge(er, on='date', how='left').merge(pm, on='date', how='left')
d['fires3'] = d.fires.rolling(3, min_periods=3).sum()                 # same day and two previous days
d['lfire'] = np.log1p(d.fires3); d['lpm'] = np.log(d.pw_pm25)
d['ym'] = d.date.dt.strftime('%Y-%m'); d['month'] = d.date.dt.month
dd = d.dropna(subset=['lpm', 'lfire', 'precip_mm', 'wind'])
fit = smf.ols('lpm ~ lfire + precip_mm + wind + C(ym)', dd).fit(cov_type='HAC', cov_kwds={'maxlags': 7})
b = fit.params['lfire']; lo, hi = fit.conf_int().loc['lfire']
dbl = lambda x: 100 * (2 ** x - 1)
mo = d.groupby('ym').agg(pm=('pw_pm25', 'mean'), fires=('fires', 'sum'), month=('month', 'first'))
son = d[d.month.isin([9, 10, 11])].dropna(subset=['pw_pm25'])
hi_fire = son.fires3 > son.fires3.median()
OUT['fires'] = {
    'days': int(len(dd)), 'b_log_fire': b, 'pct_pm_per_doubling': dbl(b), 'pct_ci': [dbl(lo), dbl(hi)], 'p': fit.pvalues['lfire'],
    'monthly_spearman': float(mo.pm.corr(mo.fires, method='spearman')),
    'son_share_of_fires': float(100 * d[d.month.isin([9, 10, 11])].fires.sum() / d.fires.sum()),
    'son_pm_high_fire_days': float(son[hi_fire].pw_pm25.mean()), 'son_pm_low_fire_days': float(son[~hi_fire].pw_pm25.mean()),
    'son_popdays_share_on_high_fire_days': float(100 * son[hi_fire].population_days_above_15.sum() / son.population_days_above_15.sum()),
}
mo.to_csv(RES / 'source_fire_monthly.csv'); d.to_csv(RES / 'source_fire_daily.csv', index=False)
print(json.dumps(OUT['fires'], indent=1))

# ---- 3. PM2.5 around the Kwekwe industrial belt -----------------------------------------------------------------
lo0, la0 = sites['Kwekwe']
rows = []
for y in range(2017, 2023):
    with rasterio.open(ROOT / f'data/public/ghap_annual_{y}.tif') as src:
        a = src.read(1).astype(float); t = src.transform
    a[(a < -9000) | (a <= 0)] = np.nan
    yy, xx = np.mgrid[0:a.shape[0], 0:a.shape[1]]
    lon = t.c + (xx + 0.5) * t.a; lat = t.f + (yy + 0.5) * t.e
    dist = np.hypot((lon - lo0) * 111.32 * np.cos(np.radians(la0)), (lat - la0) * 110.57)
    rows.append({'year': y, 'near': np.nanmean(a[dist <= 10]), 'ring': np.nanmean(a[(dist >= 20) & (dist <= 50)])})
k = pd.DataFrame(rows); k['excess'] = k.near - k.ring
OUT['kwekwe_pm25'] = {'mean_near': k.near.mean(), 'mean_ring': k.ring.mean(), 'mean_excess': k.excess.mean(),
                      'years_positive': int((k.excess > 0).sum())}
print(k.round(2).to_string(index=False))

# ---- 4. people near Manhize --------------------------------------------------------------------------------------
ee.Initialize(project=args.project)
pop = ee.ImageCollection('WorldPop/GP/100m/R2025A/v1/pop').filter(ee.Filter.calendarRange(2025, 2025, 'year')).mosaic()
pt = ee.Geometry.Point(*sites['Manhize'])
OUT['manhize_pop'] = {f'{r}km': pop.reduceRegion(ee.Reducer.sum(), pt.buffer(r * 1000), 100, maxPixels=1e9).getInfo().get('population')
                      for r in (10, 25)}
print(OUT['manhize_pop'])
json.dump(OUT, open(RES / 'source_signatures.json', 'w'), indent=1, default=float)
