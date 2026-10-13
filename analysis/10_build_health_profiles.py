"""Cumulative exposure and matched-window external health contrasts.

No local epidemiological model is fitted. RR(exp(weighted delta)) is a
representative contrast, never an aggregate rate ratio or attributable burden.
"""
from pathlib import Path
import os,json,sys
os.environ.pop('PROJ_LIB',None);os.environ.pop('PROJ_DATA',None)
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,rasterio
ROOT=Path(__file__).resolve().parents[1];RES=ROOT/'results';PUB=ROOT/'data/public';QA=ROOT/'qa'
TOTAL='Midlands Province';years=list(range(2017,2023))
plan=json.loads((ROOT/'evidence/health_profile_analysis_plan.json').read_text())
functions=pd.read_csv(ROOT/'evidence/health_profile_functions.csv')
manifest=json.loads((ROOT/'data/metadata/gee_download_manifest.json').read_text())
lookup={Path(r['file']).name:r for r in manifest}
def key(n):return n.replace(' ','_')
with np.load(QA/'spatial_support_primary.npz') as z:
    names=z['area_names'].tolist();area=z[key(TOTAL)]
with rasterio.open(PUB/'ghap_annual_2022.tif') as ds:profile=ds.profile.copy()
size=profile['height']*profile['width'];inside=area>0
with np.load(RES/'population_weights_primary_2022.npz') as z:fixed=np.stack([z[key(n)] for n in names])
with np.load(RES/'population_weights_comparator_2022.npz') as z:acfixed=np.stack([z[key(n)] for n in names])
norm=fixed/fixed.sum(axis=1)[:,None]
cum_ex=np.zeros(size);cum_days=np.zeros(size,dtype='int32');run=np.zeros(size,dtype='int32');maxrun=run.copy()
means=[];persist=np.zeros(size,dtype='int16');acmeans=[];acpersist=None
annual_pd=np.zeros(len(names));annual_expd=np.zeros(len(names));available_days=0
daily=[];lagcheck=[];previous=np.full((6,size),np.nan,dtype='float32')
for year in years:
    dates=pd.date_range(f'{year}-01-01',f'{year}-12-31');cube=np.full((len(dates),size),np.nan,dtype='float32')
    for month in range(1,13):
        file=f'ghap_daily_{year}_{month:02d}.tif';record=lookup[file]
        with rasterio.open(PUB/file) as ds:
            a=ds.read().astype('float32');a[(a==-9999)|~np.isfinite(a)]=np.nan
        ix=dates.get_indexer(pd.DatetimeIndex(record['band_dates']));cube[ix]=a.reshape(len(ix),size)
    valid=np.isfinite(cube[:,inside]).all(axis=1)
    assert np.all(np.isfinite(cube[:,inside]).any(axis=1)==valid),'Partial cell gaps need a different support policy'
    available_days+=int(valid.sum())
    excess=np.maximum(cube.astype('float64')-15,0)
    cum_ex+=np.nansum(excess,axis=0);cum_days+=(cube>15).sum(axis=0)
    for a in cube:
        run=np.where(a>15,run+1,0);maxrun=np.maximum(maxrun,run)
    am=np.nanmean(cube,axis=0,dtype='float64');means.append(am);persist+=(am>10)
    with rasterio.open(PUB/f'acag_annual_{year}.tif') as ds:
        ac=ds.read(1).astype(float).ravel();ac[ac==-9999]=np.nan
    acmeans.append(ac)
    if acpersist is None:acpersist=np.zeros(len(ac),dtype='int16')
    acpersist+=(ac>10)
    with np.load(RES/f'population_weights_primary_{year}.npz') as z:weights=np.stack([z[key(n)] for n in names])
    annual_pd+=weights@(cube>15).sum(axis=0)
    annual_expd+=weights@np.nansum(excess,axis=0)
    joined=np.concatenate([previous,cube]);ex=np.maximum(joined.astype(float)-15,0)
    # Current-date population is held fixed over each complete window, including January.
    valid_join=np.isfinite(joined[:,inside]).all(axis=1)
    wn=weights/weights.sum(axis=1)[:,None]
    exseries=np.nan_to_num(ex)@wn.T;exseries[~valid_join,:]=np.nan
    rawseries=np.nan_to_num(joined)@wn.T;rawseries[~valid_join,:]=np.nan
    fixed_ex=np.nan_to_num(ex)@norm.T;fixed_ex[~valid_join,:]=np.nan
    for j,name in enumerate(names):
        for window in [1,3,7]:
            delta=pd.Series(exseries[:,j]).rolling(window,min_periods=window).mean().to_numpy()[6:]
            raw=pd.Series(rawseries[:,j]).rolling(window,min_periods=window).mean().to_numpy()[6:]
            fx=pd.Series(fixed_ex[:,j]).rolling(window,min_periods=window).mean().to_numpy()[6:]
            daily.append(pd.DataFrame(dict(date=dates.strftime('%Y-%m-%d'),year=year,area=name,window_days=window,
                population=weights[j].sum(),mean_pm25_window=raw,delta_cap15=delta,delta_cap15_fixed2022=fx)))
        # All source reductions for the first seven dates in January validate cross-year weights independently in R.
        for day in range(7):
            for lag in range(7):
                lagcheck.append(dict(date=dates[day].strftime('%Y-%m-%d'),area=name,lag=lag,
                    excess_at_current_year_population=exseries[day+6-lag,j]))
    previous=cube[-6:].copy()
    print('Cumulative and lag profiles',year,flush=True)
    del cube,excess,joined,ex

D=pd.concat(daily,ignore_index=True);D.to_csv(RES/'health_daily_window_exposure.csv',index=False)
pd.DataFrame(lagcheck).to_csv(QA/'health_january_window_inputs.csv',index=False)
multi=np.mean(means,axis=0);acmulti=np.mean(acmeans,axis=0)
records=[]
for j,name in enumerate(names):
    w=fixed[j];aw=acfixed[j];tot=w.sum()
    rec=dict(area=name,population_2022=tot,observed_days=available_days,
        ghap_six_year_mean=float(w@np.nan_to_num(multi)/tot),acag_six_year_mean=float(aw@np.nan_to_num(acmulti)/tot),
        cumulative_excess_ug_m3_days=float(w@cum_ex/tot),mean_observed_days_above15=float(w@cum_days/tot),
        population_days_above15_annual_weights=float(annual_pd[j]),
        population_excess_concentration_days_annual_weights=float(annual_expd[j]),
        population_share_all6_above10=float(w[persist==6].sum()/tot),
        acag_population_share_all6_above10=float(aw[acpersist==6].sum()/tot),
        weighted_longest_run_days=float(w@maxrun/tot),
        maximum_cell_run_days=int(maxrun[w>0].max()))
    records.append(rec)
C=pd.DataFrame(records)
for metric in ['ghap_six_year_mean','cumulative_excess_ug_m3_days','population_days_above15_annual_weights']:
    ix=C.area!=TOTAL;C.loc[ix,'district_rank_'+metric]=C.loc[ix,metric].rank(method='min',ascending=False)
C.to_csv(RES/'cumulative_area_profiles.csv',index=False)
idx=np.flatnonzero(inside)
pd.DataFrame(dict(cell_index=idx,population_2022=fixed[names.index(TOTAL),inside],
    six_year_mean_pm25=multi[inside],cumulative_excess_ug_m3_days=cum_ex[inside],observed_days_above15=cum_days[inside],
    years_above10=persist[inside],longest_complete_run=maxrun[inside])).to_csv(RES/'cumulative_native_cell_profiles.csv',index=False)
def raster(name,a,units):
    v=a.reshape(profile['height'],profile['width']).astype('float32');v[~inside.reshape(v.shape)]=-9999
    pr=profile.copy();pr.update(count=1,dtype='float32',nodata=-9999,compress='deflate',predictor=3)
    file=RES/'rasters'/f'{name}.tif';tmp=file.with_suffix('.partial.tif')
    with rasterio.open(tmp,'w',**pr) as ds:ds.write(v,1);ds.update_tags(units=units,period='2017-2022',description='Exposure metric; not cumulative probability of disease')
    tmp.replace(file)
for name,a,unit in [('cumulative_excess_2017_2022',cum_ex,'micrograms per cubic metre days'),('years_above10_2017_2022',persist,'years'),('longest_run_2017_2022',maxrun,'days'),('mean_pm25_2017_2022',multi,'micrograms per cubic metre')]:raster(name,a,unit)

risk=[]
def emit(coef,delta,area,period,product,summary):
    row=dict(area=area,period=period,product=product,function_id=coef.id,endpoint=coef.endpoint,
        window_days=int(coef.window_days),summary=summary,delta_pm25=float(delta),reference_cap=15 if coef.window_days else 5)
    for out,col in [('rr_estimate','rr'),('ci95_lower','ci95_lower'),('ci95_upper','ci95_upper'),('pi80_lower','pi80_lower'),('pi80_upper','pi80_upper')]:
        row[out]=float(np.exp(np.log(coef[col])*delta/coef.increment_ug_m3)) if pd.notna(coef[col]) else np.nan
    risk.append(row)
for _,coef in functions.iterrows():
    if coef.window_days:
        sub=D[D.window_days==coef.window_days]
        for (name,year),g in sub.groupby(['area','year']):
            for label,q in [('median',.5),('p95',.95),('maximum',1)]:
                emit(coef,g.delta_cap15.quantile(q),name,str(year),'GHAP',label)
        for name,g in sub.groupby('area'):
            for label,q in [('median',.5),('p95',.95),('maximum',1)]:
                emit(coef,g.delta_cap15_fixed2022.quantile(q),name,'2017-2022','GHAP_fixed2022',label)
    else:
        for j,name in enumerate(names):
            for product,a,w in [('GHAP',multi,fixed[j]),('ACAG',acmulti,acfixed[j])]:
                delta=float(w@np.nan_to_num(np.maximum(a-5,0))/w.sum())
                emit(coef,delta,name,'2017-2022',product,'six_year_mean')
R=pd.DataFrame(risk);R.to_csv(RES/'outcome_specific_health_profiles.csv',index=False)
# Coefficient-aware unit tests use artificial arrays only, never study observations.
assert np.isclose((1+.021)**(10/3.28),functions.set_index('id').loc['lin2024_asthma_admissions_w7','rr_per_10'])
for w in [1,3,7]:
    a=pd.Series([10.,20.,30.,np.nan,40.,50.,60.,70.,80.,90.]);positive=(a-15).clip(lower=0)
    out=positive.rolling(w,min_periods=w).mean()
    for t in range(len(a)):
        expected=np.mean(positive.iloc[t-w+1:t+1]) if t>=w-1 and positive.iloc[t-w+1:t+1].notna().all() else np.nan
        assert (pd.isna(out[t]) and pd.isna(expected)) or np.isclose(out[t],expected)
summary=dict(period='2017-2022',available_days=available_days,functions=len(functions),risk_records=len(R),
    province=C[C.area==TOTAL].iloc[0].to_dict(),districts=C[C.area!=TOTAL].to_dict('records'),
    valid_window_dates={str(w):int(D[(D.area==TOTAL)&(D.window_days==w)].delta_cap15.notna().sum()) for w in [1,3,7]},
    acute_2022=R[(R.area==TOTAL)&(R.period=='2022')&(R.summary=='p95')].to_dict('records'),
    chronic_province=R[(R.area==TOTAL)&(R.period=='2017-2022')&(R.summary=='six_year_mean')].to_dict('records'),
    caveat='Conditional literature-transferred contrasts; no local health outcomes or event counts')
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):return [clean(v) for v in x]
    if isinstance(x,float) and not np.isfinite(x):return None
    return x
(RES/'health_profile_results.json').write_text(json.dumps(clean(summary),indent=2,allow_nan=False),encoding='utf-8')
p=RES/'analysis_summary.json';meta=json.loads(p.read_text())
meta.update(primary_scope='Population exposure, cumulative exposure and evidence-based health-risk profiles',mortality_scope='Conditional endpoint-specific main profiles and detailed supplementary scenarios',health_extension='results/health_profile_results.json')
p.write_text(json.dumps(meta,indent=2),encoding='utf-8')
print(json.dumps({k:summary[k] for k in ['available_days','functions','risk_records','valid_window_dates','province']},indent=2),flush=True)
