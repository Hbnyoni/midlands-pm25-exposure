"""Additional descriptive stability and matched-source sensitivity analyses.

Defined after the first manuscript results, at the user's request for a stronger
analysis. No local epidemiological model or data-driven source selection is used.
"""
from pathlib import Path
from datetime import datetime, timezone
import json, shutil
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];RES=ROOT/'results';E=ROOT/'evidence'
backup=ROOT/'revisions'/'before_quality_revision_2026_09_28'
backup.mkdir(parents=True,exist_ok=True)
archive=ROOT/'Frontiers_Midlands_Review_Package.zip'
if archive.exists() and not (backup/archive.name).exists():shutil.copy2(archive,backup/archive.name)
plan=dict(recorded_utc=datetime.now(timezone.utc).isoformat(),
    status='Additional analyses defined after the initial manuscript; not preregistration',
    purpose='Assess geographic stability, seasonal population-time and dependence on matched health evidence',
    rules=[
        'Retain all eight reporting areas, six years and both annual exposure products.',
        'Descending average ranks; rank 1 denotes greatest modelled exposure, not clinical risk.',
        'Compare six-year and annual concentration ranks; descriptive Spearman coefficients without inferential p-values.',
        'Retain all one-, three- and seven-day 95th-percentile daily-cap contrasts by year.',
        'Compare same-day asthma ED and admission estimates from Lin 2024 with Agache 2024 Table 2.',
        'Preserve all Agache PM2.5 lag 0-4 coefficients, including intervals spanning unity; only lag 0 transferred for a matched-source comparison.',
        'Summarise monthly population-days from the full available calendar using year-matched counts.',
        'No sum of separate lag estimates, pooled outcomes or clinical composite score.'
    ],
    source_selection=[
        dict(doi='10.1111/all.16041',url='https://onlinelibrary.wiley.com/doi/10.1111/all.16041',
             role='Matched asthma sensitivity; source Table 2 values and source-reported GRADE',
             precision='Three-decimal Table 2 values retained. ED lag 2 lower limit is 1.003 in Table 2 versus 1.004 in narrative; the table is the extraction source.'),
        dict(doi='10.17159/caj/2025/35/2.23815',url='https://www.scielo.org.za/scielo.php?lng=en&nrm=iso&pid=S2410-972X2025000200010&script=sci_arttext',
             role='Zimbabwe context only; Harare 2023-2024 measurements do not validate Midlands 2017-2022'),
        dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC8994423/',role='Screened but excluded from exact-window quantitative substitution',
             reason='COPD pooled metrics mix single-day lags or multi-day/distributed-lag windows, combine ED/admission outcomes and select strongest estimates within categories.')])
(E/'inference_revision_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')

meta=[]
for endpoint,vals in {
    'asthma_ED':[(18,1.012,1.001,1.023,'Moderate'),(15,1.013,1.004,1.022,'Moderate'),(14,1.014,1.003,1.025,'Moderate'),(14,1.023,1.007,1.039,'Moderate'),(7,1.016,.996,1.037,'Low')],
    'asthma_admissions':[(30,1.007,.992,1.023,'Low'),(31,1.006,1.000,1.012,'Low'),(31,1.011,.997,1.025,'Low'),(25,1.003,1.000,1.006,'Moderate'),(19,1.002,.997,1.007,'Low')]
}.items():
    for lag,(studies,rr,lo,hi,grade) in enumerate(vals):
        meta.append(dict(endpoint=endpoint,lag_days=lag,studies=studies,increment_ug_m3=10,
                         rr=rr,ci95_lower=lo,ci95_upper=hi,source_GRADE=grade,
                         source='Agache et al. 2024',doi='10.1111/all.16041',table='2'))
M=pd.DataFrame(meta);M.to_csv(E/'agache2024_asthma_functions.csv',index=False)
C=pd.read_csv(RES/'cumulative_area_profiles.csv');C=C[C.area!='Midlands Province'].copy()
P=pd.read_csv(RES/'exposure_product_sensitivity.csv');P=P[P.area!='Midlands Province'].copy()
A=pd.read_csv(RES/'annual_exposure_summary.csv');A=A[A.area!='Midlands Province']
D=pd.read_csv(RES/'health_daily_window_exposure.csv')
daily=pd.read_csv(RES/'daily_area_exposure_and_scenarios.csv',parse_dates=['date'])
functions=pd.read_csv(E/'health_profile_functions.csv')

annual=[]
for product,col in [('GHAP','ghap_v1_pw_pm25'),('ACAG','comparator_pw_pm25')]:
    for year,g in P.groupby('year'):
        ranks=g[col].rank(ascending=False,method='average')
        for (_,row),rank in zip(g.iterrows(),ranks):
            annual.append(dict(area=row.area,year=year,product=product,concentration=row[col],rank=rank))
AR=pd.DataFrame(annual);AR.to_csv(RES/'annual_district_rank_stability.csv',index=False)
agreement=[]
for year,g in AR.groupby('year'):
    x=g.pivot(index='area',columns='product',values='rank')
    agreement.append(dict(year=year,spearman=float(x.GHAP.corr(x.ACAG)),
                          max_absolute_rank_change=float((x.GHAP-x.ACAG).abs().max())))
G=pd.DataFrame(agreement);G.to_csv(RES/'annual_product_rank_agreement.csv',index=False)

W=D[D.area!='Midlands Province'].groupby(['area','year','window_days'],as_index=False).delta_cap15.quantile(.95)
W['rank']=W.groupby(['year','window_days']).delta_cap15.rank(ascending=False,method='average')
W.to_csv(RES/'acute_window_rank_stability.csv',index=False)
S=C[['area','ghap_six_year_mean','acag_six_year_mean','cumulative_excess_ug_m3_days','population_days_above15_annual_weights']].copy()
S['ghap_rank']=S.ghap_six_year_mean.rank(ascending=False,method='average')
S['acag_rank']=S.acag_six_year_mean.rank(ascending=False,method='average')
for label,col in [('cumulative','cumulative_excess_ug_m3_days'),('population_days','population_days_above15_annual_weights')]:
    S[label+'_rank']=S[col].rank(ascending=False,method='average')
for label,df in [('annual',AR),('acute',W)]:
    r=df.groupby('area')['rank'].agg(['min','max','median']).add_prefix(label+'_rank_')
    S=S.merge(r,on='area',validate='one_to_one')
S.to_csv(RES/'district_priority_stability.csv',index=False)

# Geographic ranks do not vary across positive log-linear health coefficients
# at the same window and contrast; evaluate magnitude and interval sensitivity.
rr=[]
for endpoint in ['asthma_ED','asthma_admissions']:
    m=M[(M.endpoint==endpoint)&(M.lag_days==0)].iloc[0]
    l=functions[(functions.endpoint==endpoint)&(functions.window_days==1)].iloc[0]
    for (area,year),g in D[D.window_days==1].groupby(['area','year']):
        delta=float(g.delta_cap15.quantile(.95))
        for source,c in [('Agache 2024',m),('Lin 2024',l)]:
            row=dict(area=area,year=year,endpoint=endpoint,source=source,window_days=1,
                     exposure_summary='p95',reference_cap=15,delta_pm25=delta)
            for output,input_ in [('rr_estimate','rr'),('ci95_lower','ci95_lower'),('ci95_upper','ci95_upper')]:
                row[output]=float(c[input_]**(delta/c.increment_ug_m3))
            rr.append(row)
RR=pd.DataFrame(rr);RR.to_csv(RES/'matched_asthma_evidence_sensitivity.csv',index=False)

daily['month']=daily.date.dt.month
T=daily.groupby(['area','month'],as_index=False).agg(
    exposed_population_days=('population_days_above_15','sum'),
    observed_population_days=('valid_exposure_population','sum'))
T['exposed_fraction']=T.exposed_population_days/T.observed_population_days
T['share_of_area_exposed_population_days']=T.exposed_population_days/T.groupby('area').exposed_population_days.transform('sum')
T.to_csv(RES/'seasonal_population_exposure.csv',index=False)
prov=T[T.area=='Midlands Province'];son=prov[prov.month.isin([9,10,11])]
summary=dict(analysis_role=plan['status'],
    six_year_spearman=float(S.ghap_rank.corr(S.acag_rank)),
    annual_spearman_min=float(G.spearman.min()),annual_spearman_max=float(G.spearman.max()),
    same_six_year_rank_count=int((S.ghap_rank==S.acag_rank).sum()),
    ghap_highest=S.loc[S.ghap_rank.idxmin(),'area'],acag_highest=S.loc[S.acag_rank.idxmin(),'area'],
    shared_six_year_upper_half=S[(S.ghap_rank<=4)&(S.acag_rank<=4)].area.tolist(),
    son_share_exposed_population_days=float(son.exposed_population_days.sum()/prov.exposed_population_days.sum()),
    son_share_observed_population_days=float(son.observed_population_days.sum()/prov.observed_population_days.sum()),
    peak_month=int(prov.loc[prov.exposed_fraction.idxmax(),'month']),
    province_2022_asthma=RR[(RR.area=='Midlands Province')&(RR.year==2022)].to_dict('records'),
    district_priorities=S.to_dict('records'),
    caveat='Descriptive input sensitivity. Model products are not independent validation; ranks are not clinical grades.')
(RES/'inference_strengthening_results.json').write_text(json.dumps(summary,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps({k:v for k,v in summary.items() if k!='district_priorities'},indent=2))
