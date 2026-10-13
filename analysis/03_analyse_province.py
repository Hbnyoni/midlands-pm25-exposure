"""Province-wide exposure assessment; mortality functions are supplementary.

Memory is bounded by reading one calendar year at a time. Calendar gaps remain
missing. The source grids and population-count conserving overlays are retained.
"""
from pathlib import Path
import json,os,sys,hashlib,platform
os.environ.pop('PROJ_LIB',None);os.environ.pop('PROJ_DATA',None)
os.environ['OPENBLAS_NUM_THREADS']='1'
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np,pandas as pd,rasterio,shapely
from risk_core import who_band,weighted_mean,transferred_rr,longest_exceedance_run
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'data/public';M=ROOT/'data/metadata';OUT=ROOT/'results';QA=ROOT/'qa';RAST=OUT/'rasters'
TOTAL='Midlands Province'
cfg=json.loads((ROOT/'analysis/study_config.json').read_text())
audit=json.loads((M/'ghap_units_audit.json').read_text())
assert audit['verified_gee_multiplier']==cfg['ghap_scale_factor']==1
manifest=json.loads((M/'gee_download_manifest.json').read_text())
lookup={Path(r['file']).name:r for r in manifest}
for year in cfg['years']:
    expected=[f'ghap_daily_{year}_{m:02d}.tif' for m in range(1,13)]+[f'{prefix}_{year}.tif' for prefix in ['ghap_annual','worldpop_population','acag_annual']]
    for name in expected:
        assert name in lookup and (P/name).exists(),f'Missing input {name}'
        with (P/name).open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==lookup[name]['sha256'],name
        if name.startswith('ghap'):assert lookup[name]['source_scale_factor_applied']==1
print('All 90 public raster hashes verified',flush=True)
def key(n):return n.replace(' ','_')
def read(path,dtype='float64'):
    with rasterio.open(path) as src:
        a=src.read().astype(dtype);a[(a==-9999)|~np.isfinite(a)]=np.nan
        if src.nodata is not None:a[a==src.nodata]=np.nan
        return a,src.profile.copy()
def same(a,b):return all(a[k]==b[k] for k in ['height','width','transform','crs'])
with np.load(QA/'spatial_support_primary.npz') as z:
    names=z['area_names'].tolist();areas={n:z[key(n)] for n in names}
with np.load(QA/'spatial_support_comparator.npz') as z:acareas={n:z[key(n)] for n in names}
_,profile=read(P/'ghap_annual_2022.tif');height,width=profile['height'],profile['width'];size=height*width
inside=areas[TOTAL]>0;indices=np.flatnonzero(inside)
rrrows,cccols=np.unravel_index(indices,(height,width));xs,ys=rasterio.transform.xy(profile['transform'],rrrows,cccols)
functions=pd.read_csv(ROOT/'evidence/published_risk_functions.csv').set_index('id')
annual_rows=[];annual_bands=[];daily_bands=[];daily_frames=[];risk_rows=[];cell_frames=[];coverage=[];products=[];consistency=[];fixed=[]
previous_day=None;missing=[];source_days=0;global_min=np.inf;global_max=-np.inf
fixed_weights={}
for year in [min(cfg['years']),max(cfg['years'])]:
    with np.load(OUT/f'population_weights_primary_{year}.npz') as z:fixed_weights[year]=z[key(TOTAL)]
def write_raster(name,values,units,description):
    a=np.asarray(values).reshape(height,width).astype('float32')
    a[~inside.reshape(height,width)]=-9999;a[~np.isfinite(a)]=-9999
    pr=profile.copy();pr.update(count=1,dtype='float32',nodata=-9999,compress='deflate',predictor=3)
    temp=RAST/(name+'.partial.tif')
    with rasterio.open(temp,'w',**pr) as dst:
        dst.write(a,1);dst.update_tags(units=units,description=description,exposure_product='GHAP V1',data_role='Public model estimate, not a local measurement')
    temp.replace(RAST/name)

for year in cfg['years']:
    dates=pd.date_range(f'{year}-01-01',f'{year}-12-31')
    yd=np.full((len(dates),size),np.nan,dtype='float32');available=np.zeros(len(dates),dtype=bool)
    for month in range(1,13):
        name=f'ghap_daily_{year}_{month:02d}.tif';record=lookup[name]
        a,pr=read(P/name,dtype='float32');assert same(profile,pr)
        dd=pd.DatetimeIndex(record['band_dates']);pos=dates.get_indexer(dd)
        assert len(dd)==a.shape[0] and np.all(pos>=0) and not available[pos].any()
        yd[pos]=a.reshape(len(pos),-1);available[pos]=True
    missing.extend(dates[~available].strftime('%Y-%m-%d').tolist());source_days+=int(available.sum())
    assert np.isfinite(yd[available][:,inside]).all(),'Spatial gaps require a revised common-support analysis'
    assert np.nanmin(yd[:,inside])>=0
    global_min=min(global_min,float(np.nanmin(yd[:,inside])));global_max=max(global_max,float(np.nanmax(yd[:,inside])))
    valid=np.isfinite(yd);fraction=valid.mean(axis=0)
    quarterly=np.stack([valid[dates.quarter==q].mean(axis=0) for q in range(1,5)])
    eligible=(fraction>=cfg['minimum_daily_coverage_fraction'])&(quarterly.min(axis=0)>=cfg['minimum_quarterly_coverage_fraction'])
    assert eligible[inside].all()
    am=np.nanmean(yd,axis=0,dtype='float64');p99=np.nanpercentile(yd,99,axis=0,method='linear')
    exceeds=yd>15;exdays=exceeds.sum(axis=0);longest=longest_exceedance_run(yd,15)
    excess=np.maximum(yd.astype('float64')-15,0)
    lag=np.full_like(excess,np.nan);lag[1:]=(excess[1:]+excess[:-1])/2
    previous_excess=None if previous_day is None else np.maximum(previous_day.astype('float64')-15,0)
    if previous_excess is not None:lag[0]=(excess[0]+previous_excess)/2
    lag_valid=np.isfinite(lag)
    annual_native,pr=read(P/f'ghap_annual_{year}.tif');assert same(profile,pr)
    diff=am-annual_native.ravel()
    consistency.append({'year':year,'maximum_absolute_daily_vs_native_annual_difference':float(abs(diff[inside]).max()),'area_mean_daily_minus_native_annual':weighted_mean(diff,areas[TOTAL])})
    with np.load(OUT/f'population_weights_primary_{year}.npz') as z:
        weights={n:z[key(n)] for n in names};centres={n:z['centre_'+key(n)] for n in names}
    with np.load(OUT/f'population_weights_comparator_{year}.npz') as z:acweights={n:z[key(n)] for n in names}
    acag,_=read(P/f'acag_annual_{year}.tif');acag=acag.ravel()
    annual_cls=who_band(am,cfg['annual_cutpoints']);day_cls=who_band(yd,cfg['daily_cutpoints'])
    fill=np.nan_to_num(yd);excess_fill=np.nan_to_num(excess);lag_fill=np.nan_to_num(lag)
    for name in names:
        w=weights[name];aw=areas[name];mask=aw>0;total=w.sum()
        valid_population=valid@w;valid_lag_population=lag_valid@w
        def series(a,den=valid_population):return np.divide(a@w,den,out=np.full(len(dates),np.nan),where=den>0)
        pw=series(fill);pd_above=exceeds@w;pwex=series(excess_fill);delta=series(lag_fill,valid_lag_population)
        prev=np.r_[np.nan if previous_excess is None else weighted_mean(previous_excess,w),pwex[:-1]]
        annual_rows.append({'year':year,'area':name,'area_km2':aw.sum()/1e6,'population':total,'population_density_per_km2':total/(aw.sum()/1e6),
            'native_cells_intersecting':int(mask.sum()),'area_weighted_pm25':weighted_mean(am,aw),'population_weighted_pm25':weighted_mean(am,w),
            'population_weighted_pm25_centre_allocation':weighted_mean(am,centres[name]),'annual_pm25_cell_min':am[mask].min(),'annual_pm25_cell_max':am[mask].max(),
            'population_above_annual_5':w[am>5].sum(),'population_above_annual_10':w[am>10].sum(),'population_above_annual_15':w[am>15].sum(),
            'population_in_cells_p99_above_15':w[p99>15].sum(),'population_fraction_p99_above_15':w[p99>15].sum()/total,
            'population_weighted_exceedance_days':weighted_mean(exdays,w),'population_days_above_15':pd_above.sum(),'total_population_days':valid_population.sum(),
            'fraction_population_days_above_15':pd_above.sum()/valid_population.sum(),'calendar_population_days':total*len(dates),
            'cell_exceedance_days_min':int(exdays[mask].min()),'cell_exceedance_days_max':int(exdays[mask].max()),
            'cell_p99_min':p99[mask].min(),'cell_p99_max':p99[mask].max(),'population_weighted_cell_p99':weighted_mean(p99,w),
            'p99_of_population_weighted_daily_series':float(np.nanpercentile(pw,99)),'area_series_days_above_15':int(np.sum(pw>15)),
            'max_fraction_population_exposed_daily_above_15':float(np.nanmax(series(exceeds))),'annual_days':len(dates),'available_days':int(np.isfinite(pw).sum()),
            'minimum_cell_completeness':fraction[mask].min(),'native_annual_product_pw_pm25':weighted_mean(annual_native.ravel(),w)})
        coverage.append({'year':year,'area':name,'calendar_days':len(dates),'minimum_daily_fraction':fraction[mask].min(),
                         'minimum_quarter_fraction':quarterly[:,mask].min(),'fraction_area_eligible':aw[eligible].sum()/aw.sum(),'fraction_population_eligible':w[eligible].sum()/total})
        for band in range(6):
            chosen=annual_cls==band
            annual_bands.append({'year':year,'area':name,'band':band,'population':w[chosen].sum(),'population_percent':100*w[chosen].sum()/total,'area_km2':aw[chosen].sum()/1e6})
            days=float(np.sum((day_cls==band)@w))
            daily_bands.append({'year':year,'area':name,'band':band,'population_days':days,'population_day_percent':100*days/valid_population.sum()})
        acute={}
        for ident,coef in functions.iterrows():
            if coef.exposure_window=='long_term':
                scenarios=[('annual',cap,weighted_mean(np.maximum(am-cap,0),w)) for cap in [5,10,15,25,35]]
            else:
                scenarios=[('lag01_'+label,15,float(np.nanquantile(delta,q))) for label,q in [('median',.5),('p95',.95),('maximum',1)]]
                for suffix,column in [('rr','rr_per_10'),('ci_lower','ci95_lower'),('ci_upper','ci95_upper')]:acute[ident+'_'+suffix]=transferred_rr(delta,coef[column])
            for window,cap,change in scenarios:
                values={k:float(transferred_rr(change,coef[k])) if pd.notna(coef[k]) else np.nan for k in ['rr_per_10','ci95_lower','ci95_upper','pi80_lower','pi80_upper']}
                values['rr_estimate']=values.pop('rr_per_10')
                risk_rows.append({'year':year,'area':name,'function_id':ident,'endpoint':coef.outcome,'window':window,'reference_cap':cap,'delta_pm25':change,**values})
        daily_frames.append(pd.DataFrame({'date':dates.strftime('%Y-%m-%d'),'year':year,'area':name,'population':total,
            'valid_exposure_population':valid_population,'pw_pm25':pw,'population_fraction_above_15':series(exceeds),'population_days_above_15':pd_above,
            'pw_daily_excess_above_15':pwex,'pw_previous_day_excess_at_current_population':prev,'pw_lag01_reduction':delta,**acute}))
        acw=acweights[name];acmask=acareas[name]>0
        assert np.isfinite(acag[acmask]).all() and abs(acw.sum()-total)<.01
        acmean=weighted_mean(acag,acw);ghmean=weighted_mean(am,w)
        products.append({'year':year,'area':name,'product':'ACAG V6.GL.02','ghap_v1_pw_pm25':ghmean,'comparator_pw_pm25':acmean,
            'difference_comparator_minus_ghap':acmean-ghmean,'comparator_area_mean':weighted_mean(acag,acareas[name]),
            'comparator_population':acw.sum(),'comparator_population_above_5':acw[acag>5].sum(),'comparator_population_above_10':acw[acag>10].sum(),'comparator_population_above_15':acw[acag>15].sum()})
    chronic=transferred_rr(np.maximum(am-5,0),functions.loc['orellano2024_respiratory','rr_per_10'])
    cell_frames.append(pd.DataFrame({'year':year,'cell_index':indices,'row':rrrows,'col':cccols,'longitude':xs,'latitude':ys,
        'clipped_area_km2':areas[TOTAL][inside]/1e6,'population':weights[TOTAL][inside],'population_centre':centres[TOTAL][inside],
        'annual_pm25':am[inside],'annual_who_band':annual_cls[inside],'daily_p99':p99[inside],'days_above_15':exdays[inside],
        'longest_exceedance_run':longest[inside],'annual_rr_respiratory_scenario':chronic[inside]}))
    for file,value,unit,description in [
        ('annual_pm25',am,'micrograms/m3','Available-day annual mean'),('annual_who_band',annual_cls,'WHO interval 0 to 5','Upper-inclusive WHO annual exposure interval'),
        ('daily_p99',p99,'micrograms/m3','Annual daily 99th percentile'),('days_above_15',exdays,'days','Observed daily estimates strictly above 15 micrograms/m3'),
        ('longest_run_above_15',longest,'days','Longest within-year run of observed exceedances'),('population_on_pm_grid',weights[TOTAL],'persons per clipped exposure cell','Count-conserving native WorldPop allocation'),
        ('chronic_respiratory_rr_scenario',chronic,'relative contrast','Supplementary literature-transferred contrast; not local mortality')]:
        write_raster(f'{file}_{year}.tif',value,unit,description)
    for py,w in fixed_weights.items():fixed.append({'exposure_year':year,'fixed_population_year':py,'pw_pm25':weighted_mean(am,w)})
    previous_day=yd[-1].copy()
    print('Analysed',year,'population',round(weights[TOTAL].sum()),'PM2.5',round(weighted_mean(am,weights[TOTAL]),4),flush=True)
    del yd,valid,excess,lag,lag_valid,fill,excess_fill,lag_fill,day_cls,exceeds

for records,name in [(annual_rows,'annual_exposure_summary'),(annual_bands,'annual_who_population_bands'),(daily_bands,'daily_who_population_day_bands'),
    (risk_rows,'transferred_rr_scenarios'),(coverage,'exposure_completeness'),(products,'exposure_product_sensitivity'),(fixed,'fixed_population_sensitivity')]:
    pd.DataFrame(records).to_csv(OUT/(name+'.csv'),index=False)
daily=pd.concat(daily_frames,ignore_index=True);daily.to_csv(OUT/'daily_area_exposure_and_scenarios.csv',index=False)
pd.concat(cell_frames,ignore_index=True).to_csv(OUT/'native_cell_annual_metrics.csv',index=False)
pd.DataFrame(consistency).to_csv(QA/'ghap_daily_annual_consistency.csv',index=False)
daily['month']=pd.to_datetime(daily.date).dt.month
daily.groupby(['area','year','month'],as_index=False).agg(pm25_mean=('pw_pm25','mean'),pm25_max=('pw_pm25','max'),
    mean_population_fraction_above_15=('population_fraction_above_15','mean'),calendar_days=('date','size'),available_days=('pw_pm25','count')).to_csv(OUT/'monthly_exposure.csv',index=False)
summary={'calendar_start':f'{min(cfg["years"])}-01-01','calendar_end':f'{max(cfg["years"])}-12-31','calendar_days':source_days+len(missing),
    'available_daily_grids':source_days,'missing_daily_dates':missing,'source_grid_crs':str(profile['crs']),'source_grid_transform':list(profile['transform'])[:6],
    'native_grid_shape':[height,width],'native_cells_intersecting_province':int(inside.sum()),'study_area_km2':float(areas[TOTAL].sum()/1e6),
    'minimum_daily_pm25':global_min,'maximum_daily_pm25':global_max,'reporting_areas':names,'primary_scope':'Population exposure and product sensitivity',
    'mortality_scope':'Supplementary transferred scenarios only','local_health_data_used':False,'synthetic_exposure_used':False,
    'missing_policy':'No imputation; source dates retained; missing days excluded from valid denominators and retained in lag calendar',
    'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'rasterio':rasterio.__version__,'shapely':shapely.__version__,
    'primary_annual_estimate':'Mean of available daily grids after 90% annual and 75% quarterly completeness checks',
    'uncertainty_note':'Product and allocation sensitivities; published risk intervals exclude local exposure and population model uncertainty'}
(OUT/'analysis_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2),flush=True)
