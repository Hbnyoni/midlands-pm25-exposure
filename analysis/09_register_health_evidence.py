"""Record endpoint and window choices before the additional spatial calculations."""
from pathlib import Path
import json,csv,math
ROOT=Path(__file__).resolve().parents[1]
rows=[]
def add(ident,source,endpoint,window,increment,rr,lo,hi,pi_lo=None,pi_hi=None,population='',notes=''):
    rows.append(dict(id=ident,source=source,endpoint=endpoint,window_days=window,
        increment_ug_m3=increment,rr=rr,ci95_lower=lo,ci95_upper=hi,
        pi80_lower=pi_lo,pi80_upper=pi_hi,rr_per_10=rr**(10/increment),
        ci95_lower_per_10=lo**(10/increment),ci95_upper_per_10=hi**(10/increment),
        population=population,notes=notes))
for endpoint,values in {
    'asthma_admissions':[(1,5.10,1.7,.8,2.7),(3,4.09,2.1,1.2,3.0),(7,3.28,2.1,1.1,3.2)],
    'asthma_ED':[(1,4.88,2.1,1.7,2.4),(3,3.96,1.9,1.5,2.2),(7,3.05,2.4,2.0,2.8)],
    'COPD_admissions':[(1,4.99,1.5,.5,2.4),(3,3.97,1.9,1.0,2.8),(7,3.11,3.1,2.0,4.2)],
    'COPD_ED':[(1,4.84,1.2,.2,2.4),(3,3.91,1.3,.2,2.3),(7,3.06,2.4,1.2,3.6)]
}.items():
    for w,iqr,er,lo,hi in values:
        add(f'lin2024_{endpoint}_w{w}','10.1016/j.envpol.2024.124585',endpoint,w,iqr,
            1+er/100,1+lo/100,1+hi/100,population='Residents of all ages within 15 miles of six New York State monitors, 2014-2019',
            notes='Table 1 total PM2.5 mass; full-period estimate. Each window is a separate moving-average model. Source interprets conditional logistic estimate as a rate ratio. No component-specific coefficients used.')
for endpoint,er,lo,hi in [('cardiovascular_admissions',.26,.17,.35),('IHD_admissions',.31,.22,.40),('ischaemic_stroke_admissions',.29,.18,.40)]:
    add(f'tian2019_{endpoint}_w1','10.1136/bmj.l6572',endpoint,1,10,1+er/100,1+lo/100,1+hi/100,
        population='Adults covered by urban employee insurance in 184 Chinese cities, 2014-2017',notes='Same-day total PM2.5; pooled log-linear approximation. Cause-specific endpoints overlap the cardiovascular total.')
add('tian2019_pneumonia_admissions_w3','10.1371/journal.pmed.1003010','pneumonia_admissions',3,10,1.0031,1.0015,1.0046,
    population='Adults covered by urban employee insurance in 184 Chinese cities, 2014-2017',notes='Three-day mean lag 0-2. Primary single-pollutant estimate; NO2-adjusted estimate attenuated and CI crossed unity.')
for outcome,rr,lo,hi,pl,pu in [
    ('all_cause_mortality',1.095,1.064,1.127,.966,1.241),
    ('circulatory_mortality',1.127,1.102,1.152,1.046,1.214),
    ('IHD_mortality',1.143,1.102,1.186,1.026,1.273),
    ('cerebrovascular_mortality',1.146,1.101,1.192,1.034,1.270),
    ('ALRI_mortality',1.204,1.095,1.325,1.005,1.443),
    ('lung_cancer_mortality',1.093,1.053,1.135,.999,1.196),
    ('respiratory_mortality',1.136,1.079,1.197,.990,1.305),
    ('COPD_mortality',1.138,1.080,1.198,1.023,1.266)]:
    add('orellano2024_'+outcome,'10.3389/ijph.2024.1607683',outcome,0,10,rr,lo,hi,pl,pu,
        population='Endpoint-specific cohort populations in WHO guideline evidence update',
        notes='Table 1 PM2.5 mortality; long-term association, not disease incidence or a six-year cumulative probability. Prediction interval is 80%.')
out=ROOT/'evidence/health_profile_functions.csv'
with out.open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
plan={
 'recorded_date':'2026-09-24','reason':'User requested health-risk profiling as a core analysis, including published health endpoints and lag windows.',
 'selection_type':'Targeted primary-source evidence selection, not a systematic review or complete literature census',
 'selection_rule':'Total ambient PM2.5, exact reported endpoint, exposure increment and averaging window, retrievable confidence limits; long-term meta-analysis plus multi-city morbidity studies.',
 'acute_windows':'Lin 2024 same-day, current plus previous two days, current plus previous six days; Tian cardiovascular same-day; Tian pneumonia current plus previous two days. All Lin windows 1, 3 and 7 selected together before additional spatial calculations, not according to local significance or maximum effect.',
 'reference_scenario':'Cap each daily cell concentration at 15 before calculating complete moving-average reductions. Annual or six-year mean cell concentrations capped at 5. References are WHO comparison levels, not biological thresholds.',
 'cumulative_exposure':'Observed-day positive excess concentration above 15 summed across 2017-2022 in micrograms per cubic metre days. Fixed 2022 population for spatial intensity comparisons; year-matched population for historical population-days. Missing days break runs and invalidate every window containing them.',
 'chronic_profile':'Six-year mean of eligible annual concentrations with fixed 2022 population on each product native grid. Representative exposure contrast evaluated with each endpoint function; not average person-level risk or disease probability.',
 'classification':'WHO concentration intervals, number of years above annual interim target 4 (10), continuous cumulative excess and descriptive district ranks. No arbitrary clinically labelled composite risk score.',
 'excluded_or_context_only':[
  {'source':'10.1136/bmj-2023-076322','reason':'Sun et al. 2024 headline coefficients are restricted to PM2.5 below 15; not transferred to daily excess above 15. Supports context that guideline level is not a no-effect boundary.'},
  {'source':'10.1016/j.envpol.2024.124585','reason':'POC, SOC and source-specific coefficients excluded because GHAP estimates total ambient PM2.5 mass.'},
  {'source':'10.1136/bmj.l6572','reason':'Haemorrhagic stroke estimate was -0.02% (95% CI -0.23 to 0.19) per 10 and did not support a positive association. No general claim that every stroke subtype increases.'}
 ],
 'interpretation_limits':['No local health outcomes, baseline rates, individual risks or attributable event counts','All-age residential weights locate exposure; they do not make adult or disease-specific coefficients applicable to all residents','Do not add overlapping endpoints or multiply daily ratios into cumulative disease risk','Separate moving-average fits do not identify independent day-specific lag contributions','Coefficient intervals exclude environmental model and population uncertainty; product sensitivity is separate']}
(ROOT/'evidence/health_profile_analysis_plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
cfgpath=ROOT/'analysis/study_config.json';cfg=json.loads(cfgpath.read_text())
cfg.update(design='Historical ecological exposure assessment and evidence-based health-risk profiling',
 primary_transferred_endpoint='Outcome-specific conditional morbidity and mortality contrasts',
 secondary_transferred_endpoint='Product and averaging-window sensitivities',
 mortality_analysis_role='Long-term outcome-specific contrasts in main profiles; detailed annual and acute mortality scenarios supplementary',
 health_profile_plan='evidence/health_profile_analysis_plan.json',acute_profile_window_lengths=[1,3,7])
cfg['primary_results'] = list(dict.fromkeys(cfg['primary_results']+['Cumulative and persistent exposure profiles','Published outcome-specific health contrasts with matched exposure windows']))
cfgpath.write_text(json.dumps(cfg,indent=2),encoding='utf-8')
print(f'Recorded {len(rows)} matched functions before health-profile calculations')
