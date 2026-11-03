"""Source-attributed exposure model for Midlands Province, 2017 to 2022 (3 November 2026).

Every cell-day of GHAP PM2.5 (C) is split into three additive parts:
  fire-attributable   C - C0,  C0 = C (1 + F)^-beta, where F is MODIS fire pixels within 50 km over the day and the two before
  industrial/urban    I, the local excess of the cell's annual mean over a 20 to 50 km ring, within 10 km of a NO2 source
                      whose 2018 to 2024 NO2 excess is above zero (95% interval), applied to every day of that year
  background          C0 - I (not below zero)
beta comes from a district-day panel, ln C_jd = beta ln(1 + F_jd) + g1 R_d + g2 U_d + a_(j, year-month), with standard errors
clustered by year-month. Exceedance population-days (C > 15) are attributed with the Shapley average of the two removal orders.
A monthly version of the fire model is fitted to GHAP and to ACAG V5.GL monthly PM2.5 so the fire effect is estimated under
both products. Outputs in results/model_*.
Usage: python analysis/34_source_model.py   (after 28, 29 and 33)
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from scipy import ndimage
import statsmodels.formula.api as smf
from shapely.geometry import shape

ROOT = Path(__file__).resolve().parents[1]
PUB, SRC, RES = ROOT / 'data/public', ROOT / 'data/sources', ROOT / 'results'
YEARS = range(2017, 2023)
DIST = ['Chirumhanzu', 'Gokwe North', 'Gokwe South', 'Gweru', 'Kwekwe', 'Mberengwa', 'Shurugwi', 'Zvishavane']
key = lambda n: n.replace(' ', '_')

# ---------------------------------------------------------------- grids
with rasterio.open(PUB / 'ghap_daily_2017_01.tif') as s:
    T, SHP = s.transform, s.shape
rows, cols = np.mgrid[0:SHP[0], 0:SHP[1]]
LON, LAT = T.c + (cols + .5) * T.a, T.f + (rows + .5) * T.e
FW, FN, FS = 27.0, -16.4, 0.05
PIX = 28.6   # the grid stores the burning fraction of each 0.05 degree cell; x 28.6 = 1 km fire pixels (checked against FIRMS totals)
fi = np.clip(((FN - LAT) / FS).astype(int), 0, 109); fj = np.clip(((LON - FW) / FS).astype(int), 0, 109)
def disc(radius_km):
    ry, rx = radius_km / (FS * 110.57), radius_km / (FS * 111.32 * np.cos(np.radians(19)))
    yy, xx = np.mgrid[-int(ry) - 1:int(ry) + 2, -int(rx) - 1:int(rx) + 2]
    return ((yy / ry) ** 2 + (xx / rx) ** 2 <= 1).astype(float)
KER = {r: disc(r) for r in (25, 50, 100)}
W = {y: np.load(RES / f'population_weights_primary_{y}.npz') for y in YEARS}
POP = {y: W[y]['Midlands_Province'].reshape(SHP) for y in YEARS}
feats = json.loads((PUB / 'study_administrative_units.geojson').read_text())['features']
AREA = {f['properties']['shapeName']: ~geometry_mask([shape(f['geometry'])], SHP, T) for f in feats}

def fire_days():
    """Yield (date, {radius: day grid}) for every day 2016-12-30 onward, fire counts convolved with each radius."""
    for y in YEARS:
        for m in range(1, 13):
            with rasterio.open(SRC / f'fire_grid_{y}_{m:02d}.tif') as s:
                a = s.read().astype(np.float32) * PIX; names = s.descriptions
            for k in range(a.shape[0]):
                date = pd.Timestamp(names[k].split('_')[-1] if names[k] else f'{y}-{m:02d}-{k + 1:02d}')
                yield date, {r: ndimage.convolve(a[k], K, mode='constant') for r, K in KER.items()}

def ghap_days():
    for y in YEARS:
        for m in range(1, 13):
            with rasterio.open(PUB / f'ghap_daily_{y}_{m:02d}.tif') as s:
                a = s.read().astype(np.float32); names = s.descriptions
            a[(a <= -9000) | ~np.isfinite(a)] = np.nan
            for k in range(a.shape[0]): yield pd.Timestamp(names[k]), a[k]

# ---------------------------------------------------------------- pass 1: district daily series and fire covariates
fires = {d: g for d, g in fire_days()}
dates = sorted(fires)
def window(date, lag, r):
    """Fire pixels within r km summed over `lag` days ending on date, on the fire grid."""
    out = np.zeros_like(fires[dates[0]][r])
    for k in range(lag):
        g = fires.get(date - pd.Timedelta(days=k))
        if g is not None: out = out + g[r]
    return out
def lead(date, r):
    out = np.zeros_like(fires[dates[0]][r])
    for k in (3, 4, 5):
        g = fires.get(date + pd.Timedelta(days=k))
        if g is not None: out = out + g[r]
    return out

CACHE = RES / 'model_district_daily.csv'
rows_ = []
annual_sum = {y: np.zeros(SHP) for y in YEARS}; annual_n = {y: np.zeros(SHP) for y in YEARS}
for date, C in ([] if CACHE.exists() else ghap_days()):
    y = date.year
    wins = {f'f{r}_{lag}': window(date, lag, r)[fi, fj] for r in (25, 50, 100) for lag in (1, 3, 7)}
    wins['lead50'] = lead(date, 50)[fi, fj]
    for d in DIST:
        w = W[y][key(d)].reshape(SHP); ok = np.isfinite(C) & (w > 0)
        if not ok.any(): continue
        row = {'date': date, 'area': d, 'pm': float((w[ok] * C[ok]).sum() / w[ok].sum()),
               'pm_area': float(np.nanmean(C[AREA[d]]))}
        for name, f in wins.items(): row[name] = float((w[ok] * f[ok]).sum() / w[ok].sum())
        rows_.append(row)
    annual_sum[y] += np.nan_to_num(C); annual_n[y] += np.isfinite(C)
if CACHE.exists():
    panel = pd.read_csv(CACHE, parse_dates=['date']).drop(columns=['precip_mm', 'wind', 'ym', 'g', 'lpm', 'month'])
    z = np.load(RES / 'model_annual_ghap.npz'); ANN = {y: z[str(y)] for y in YEARS}
else:
    panel = pd.DataFrame(rows_)
    ANN = {y: annual_sum[y] / np.where(annual_n[y] > 0, annual_n[y], np.nan) for y in YEARS}
    np.savez_compressed(RES / 'model_annual_ghap.npz', **{str(y): ANN[y] for y in YEARS})
era = pd.read_csv(SRC / 'era5_daily.csv', parse_dates=['date'])
panel = panel.merge(era, on='date', how='left')
panel['ym'] = panel.date.dt.strftime('%Y-%m'); panel['g'] = panel.area + '_' + panel.ym
panel['lpm'] = np.log(panel.pm); panel['month'] = panel.date.dt.month
panel.to_csv(CACHE, index=False)
print('panel', panel.shape)

# ---------------------------------------------------------------- fire model and sensitivity
def within(df, cols_, grp='g'):
    return df[cols_] - df.groupby(grp)[cols_].transform('mean')
def fit(df, xvar, extra=(), grp='g', cluster='ym'):
    d = df.dropna(subset=['lpm', xvar, 'precip_mm', 'wind', *extra]).copy()
    d['lx'] = np.log1p(d[xvar]); xs = ['lx', 'precip_mm', 'wind'] + [f'l_{e}' for e in extra]
    for e in extra: d[f'l_{e}'] = np.log1p(d[e])
    Z = within(d, ['lpm'] + xs, grp)
    y_, X = Z['lpm'].to_numpy(), Z[xs].to_numpy()
    XtX = np.linalg.inv(X.T @ X); b = XtX @ X.T @ y_; u = y_ - X @ b
    meat = np.zeros((len(xs), len(xs)))
    for _, idx in d.groupby(cluster).indices.items():
        s = X[idx].T @ u[idx]; meat += np.outer(s, s)
    G = d[cluster].nunique(); n, k = X.shape
    V = XtX @ meat @ XtX * G / (G - 1) * (n - 1) / (n - k - d[grp].nunique())
    se = np.sqrt(np.diag(V))
    return {x: (b[i], se[i]) for i, x in enumerate(xs)}, len(d)
pct = lambda b: 100 * (2 ** b - 1)
coef = []
def add(label, res, n, name='lx'):
    b, se = res[name]; coef.append({'model': label, 'beta': b, 'se': se, 'pct_per_doubling': pct(b),
                                     'lo': pct(b - 1.96 * se), 'hi': pct(b + 1.96 * se), 'n': n})
res, n = fit(panel, 'f50_3'); add('Main: 50 km, 3 days', res, n); BETA = res['lx'][0]
for r in (25, 100): res_, n_ = fit(panel, f'f{r}_3'); add(f'Radius {r} km', res_, n_)
for lag in (1, 7): res_, n_ = fit(panel, f'f50_{lag}'); add(f'{lag}-day window', res_, n_)
res_, n_ = fit(panel, 'f50_3', extra=('lead50',)); add('Placebo: fires 3 to 5 days later', res_, n_, name='l_lead50'); add('Main, with later fires held fixed', res_, n_)
son = panel[panel.month.isin([9, 10, 11])]; res_, n_ = fit(son, 'f50_3'); add('September to November only', res_, n_)
oth = panel[~panel.month.isin([9, 10, 11])]; res_, n_ = fit(oth, 'f50_3'); add('Other months', res_, n_)
for d in DIST:
    res_, n_ = fit(panel[panel.area == d], 'f50_3'); add(f'District: {d}', res_, n_)
coef = pd.DataFrame(coef); coef.to_csv(RES / 'model_fire_coefficients.csv', index=False)
print(coef.round(3).to_string(index=False))

# monthly model under both products
acag = pd.read_csv(SRC / 'acag_monthly_districts.csv')
mon = panel.groupby(['area', 'ym']).agg(pm_area=('pm_area', 'mean'), fire=('f50_1', 'sum'), precip_mm=('precip_mm', 'sum'),
                                        wind=('wind', 'mean')).reset_index()
mon = mon.merge(acag.rename(columns={'month': 'ym', 'pm25': 'acag'}), on=['area', 'ym'], how='left')
mon['year'] = mon.ym.str[:4]; mon['cal'] = mon.ym.str[5:]; mon['jc'] = mon.area + '_' + mon.cal
mono = []
for prod, col in (('GHAP', 'pm_area'), ('ACAG', 'acag')):
    m_ = mon.dropna(subset=[col]).copy(); m_['lpm'] = np.log(m_[col]); m_['lf'] = np.log1p(m_.fire)
    # string columns enter as fixed effects
    f_ = smf.ols('lpm ~ lf + precip_mm + wind + jc + year', m_).fit(cov_type='cluster', cov_kwds={'groups': m_.ym})
    b, se = f_.params['lf'], f_.bse['lf']
    mono.append({'product': prod, 'beta': b, 'se': se, 'pct_per_doubling': pct(b), 'lo': pct(b - 1.96 * se), 'hi': pct(b + 1.96 * se), 'n': len(m_)})
mono = pd.DataFrame(mono); mono.to_csv(RES / 'model_fire_monthly_products.csv', index=False); print(mono.round(3).to_string(index=False))

# ---------------------------------------------------------------- industrial and urban local excess fields
sites = json.loads((SRC / 'hotspots.json').read_text())
ms = pd.read_csv(SRC / 'monthly_sites.csv'); ms['ex'] = ms.no2_near - ms.no2_ring
keep = []
for s_, g in ms[ms.month < '2024-06'].groupby('site'):
    g = g.dropna(subset=['ex']).reset_index(drop=True)
    f_ = smf.ols('ex ~ 1', g).fit(cov_type='HAC', cov_kwds={'maxlags': 6})
    lo = f_.conf_int().iloc[0, 0]; keep.append({'site': s_, 'no2_excess': f_.params.iloc[0], 'lo': lo, 'source': bool(lo > 0)})
keep = pd.DataFrame(keep); keep.to_csv(RES / 'model_no2_sources.csv', index=False); print(keep.round(3).to_string(index=False))
SOURCES = [s_ for s_ in keep[keep.source].site if s_ != 'Manhize']   # Manhize lies outside the province
def km(lo0, la0): return np.hypot((LON - lo0) * 111.32 * np.cos(np.radians(la0)), (LAT - la0) * 110.57)
near = {s_: km(*sites[s_]) <= 10 for s_ in SOURCES}
any_near = np.zeros(SHP, bool)
for m_ in near.values(): any_near |= m_
def excess_field(A):
    I = np.zeros(SHP); rows2 = []
    for s_ in SOURCES:
        dd = km(*sites[s_]); ring = (dd >= 20) & (dd <= 50) & ~any_near & np.isfinite(A)
        base = np.nanmean(A[ring]); inc = np.where(near[s_] & np.isfinite(A), np.clip(A - base, 0, None), 0)
        I = np.maximum(I, inc); rows2.append((s_, np.nanmean(A[near[s_]]) - base))
    return I, rows2
IND, ind_rows = {}, []
acag_ann = {}
for y in YEARS:
    IND[y], r1 = excess_field(ANN[y])
    with rasterio.open(PUB / f'acag_annual_{y}.tif') as s:
        a = s.read(1).astype(float); a[(a <= -9000) | (a <= 0)] = np.nan
    acag_ann[y] = np.full(SHP, np.nan); acag_ann[y][:a.shape[0], :a.shape[1]] = a
    _, r2 = excess_field(acag_ann[y])
    ind_rows += [{'year': y, 'site': s_, 'ghap_excess': e1, 'acag_excess': e2} for (s_, e1), (_, e2) in zip(r1, r2)]
pd.DataFrame(ind_rows).to_csv(RES / 'model_industrial_excess.csv', index=False)
print(pd.DataFrame(ind_rows).groupby('site')[['ghap_excess', 'acag_excess']].agg(['mean', 'min', 'max']).round(2))

# ---------------------------------------------------------------- pass 2: decomposition and attribution
acc = {k: {d: {y: 0.0 for y in YEARS} for d in DIST + ['Midlands Province']} for k in
       ('pd_total', 'pd_fire', 'pd_ind', 'pd_bg', 'pd_obs', 'c', 'cf', 'ci', 'cb', 'wsum')}
cell = {k: np.zeros(SHP) for k in ('ex_total', 'ex_fire', 'fire_inc_son', 'n_son', 'n')}
daily = []
for date, C in ghap_days():
    y = date.year; F = window(date, 3, 50)[fi, fj]
    C0 = C * (1 + F) ** (-BETA); I = IND[y]; B = np.clip(C0 - I, 0, None)
    e_t, e_f, e_b = C > 15, C0 > 15, B > 15
    e_i_only = (C - np.minimum(I, C)) > 15                     # industry removed first
    fire_sh = .5 * ((e_t & ~e_f) + (e_i_only & ~e_b)).astype(float)   # Shapley: fire removed first, then second
    ind_sh = .5 * ((e_f & ~e_b) + (e_t & ~e_i_only)).astype(float)
    ok = np.isfinite(C)
    for d in DIST + ['Midlands Province']:
        w = W[y][key(d)].reshape(SHP) * ok
        a = acc
        a['pd_obs'][d][y] += w.sum(); a['pd_total'][d][y] += (w * e_t).sum(); a['pd_fire'][d][y] += (w * fire_sh).sum()
        a['pd_ind'][d][y] += (w * ind_sh).sum(); a['pd_bg'][d][y] += (w * e_b).sum()
        for k_, arr in (('c', C), ('cf', C - C0), ('ci', C0 - B), ('cb', B)): a[k_][d][y] += np.nansum(w * arr)
        a['wsum'][d][y] += w.sum()
    cell['ex_total'] += np.nan_to_num(e_t); cell['ex_fire'] += np.nan_to_num(fire_sh); cell['n'] += ok
    if date.month in (9, 10, 11): cell['fire_inc_son'] += np.nan_to_num(C - C0); cell['n_son'] += ok
    w = W[y]['Midlands_Province'].reshape(SHP) * ok
    daily.append({'date': date, 'c': np.nansum(w * C) / w.sum(), 'fire': np.nansum(w * (C - C0)) / w.sum(),
                  'ind': np.nansum(w * (C0 - B)) / w.sum(), 'bg': np.nansum(w * B) / w.sum()})
pd.DataFrame(daily).to_csv(RES / 'model_daily_province.csv', index=False)
out = []
for d in DIST + ['Midlands Province']:
    for y in YEARS:
        ws = acc['wsum'][d][y]
        out.append({'area': d, 'year': y, 'pm': acc['c'][d][y] / ws, 'fire': acc['cf'][d][y] / ws, 'industry': acc['ci'][d][y] / ws,
                    'background': acc['cb'][d][y] / ws, 'pd_obs': acc['pd_obs'][d][y], 'pd_total': acc['pd_total'][d][y],
                    'pd_fire': acc['pd_fire'][d][y], 'pd_industry': acc['pd_ind'][d][y],
                    'pd_other': acc['pd_total'][d][y] - acc['pd_fire'][d][y] - acc['pd_ind'][d][y]})
dy = pd.DataFrame(out); dy.to_csv(RES / 'model_district_year.csv', index=False)
np.savez_compressed(RES / 'model_cells.npz', fire_share=np.where(cell['ex_total'] > 0, cell['ex_fire'] / cell['ex_total'], np.nan),
                    fire_inc_son=cell['fire_inc_son'] / np.where(cell['n_son'] > 0, cell['n_son'], np.nan),
                    ex_days=cell['ex_total'] / 6, fire_days=cell['ex_fire'] / 6,
                    ind_mean=np.mean([IND[y] for y in YEARS], 0), pm_mean=np.nanmean([ANN[y] for y in YEARS], 0),
                    acag_mean=np.nanmean([acag_ann[y] for y in YEARS], 0), pop2022=POP[2022])
tot = dy.groupby('area')[['pd_total', 'pd_fire', 'pd_industry', 'pd_other', 'pd_obs']].sum()
tot['fire_share'] = 100 * tot.pd_fire / tot.pd_total; tot['ind_share'] = 100 * tot.pd_industry / tot.pd_total
means = dy.groupby('area')[['pm', 'fire', 'industry', 'background']].mean()
summary = tot.join(means); summary.to_csv(RES / 'model_district_summary.csv')
print(summary.round(2).to_string())
json.dump({'beta': BETA, 'sources': SOURCES}, open(RES / 'model_meta.json', 'w'), indent=1)
