"""Final figure set: readable at journal width, with explicit sensitivity.

Run after 05_make_figures.py and 15_strengthen_inference.py. Panel letters are
the only titles; interpretation and study windows belong in DOCX captions.
"""
from figure_style import *
from matplotlib.lines import Line2D
import calendar

C=pd.read_csv(RES/'cumulative_area_profiles.csv')
R=pd.read_csv(RES/'outcome_specific_health_profiles.csv')
H=pd.read_csv(RES/'matched_asthma_evidence_sensitivity.csv')
S=pd.read_csv(RES/'district_priority_stability.csv')
F=pd.read_csv(ROOT/'evidence/health_profile_functions.csv')
M=pd.read_csv(ROOT/'evidence/agache2024_asthma_functions.csv')
order=C[C.area!='Midlands Province'].sort_values('cumulative_excess_ug_m3_days',ascending=False).area.tolist()
product_handles=[Line2D([],[],marker='o',ls='',color=teal,label='GHAP'),Line2D([],[],marker='s',ls='',color=orange,label='ACAG')]

def panel(ax,letter):
    ax.set_title(letter,loc='left',fontweight='bold',pad=9)
    ax.tick_params(labelsize=8)
    ax.grid(axis='x',color='#dce3e5',lw=.6,alpha=.7);ax.set_axisbelow(True)

def contrast_text(im,value):
    rgb=np.asarray(im.cmap(im.norm(value))[:3]);linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
    return 'white' if float(linear@np.array([.2126,.7152,.0722]))<.22 else navy

def horizontal_key(fig,im,pos,label,ticks=None):
    cb=fig.colorbar(im,cax=fig.add_axes(pos),orientation='horizontal',ticks=ticks)
    cb.set_label(label,fontsize=9);cb.ax.tick_params(labelsize=8,length=2)
    cb.outline.set_linewidth(.5)

# Daily series plus monthly population-time, rather than a second concentration summary.
fig=plt.figure(figsize=(7.2,5.5));ax=fig.add_axes([.11,.57,.85,.34]);panel(ax,'a')
q=daily[daily.area=='Midlands Province']
ax.plot(q.date,q.pw_pm25,lw=.65,color=teal)
ax.axhline(15,color=orange,lw=1,ls=(0,(4,3)))
ax.set_ylabel('PM₂.₅ (µg/m³)',fontsize=9)
ax.set_xlim(pd.Timestamp('2017-01-01'),pd.Timestamp('2022-12-31'))
ax.grid(axis='x',visible=False);ax.grid(axis='y',color='#e3e8ea',lw=.6)
fig.legend(handles=[Line2D([],[],color=teal,lw=1,label='Daily population-weighted mean'),Line2D([],[],color=orange,ls='--',label='15 µg/m³')],
           loc='upper center',bbox_to_anchor=(.54,.995),ncol=2,frameon=False,fontsize=8)
ax=fig.add_axes([.11,.17,.75,.27]);panel(ax,'b');ax.grid(False)
m=monthly[monthly.area=='Midlands Province'].pivot(index='year',columns='month',values='mean_population_fraction_above_15')*100
im=ax.imshow(m,aspect='auto',cmap='YlGnBu',vmin=0,vmax=100)
ax.set_xticks(range(12),[calendar.month_abbr[i] for i in range(1,13)],fontsize=8)
ax.set_yticks(range(6),m.index,fontsize=8);ax.tick_params(length=0,pad=6)
for (i,j),v in np.ndenumerate(m.to_numpy()):ax.text(j,i,f'{v:.0f}',ha='center',va='center',fontsize=8,color=contrast_text(im,v))
cb=fig.colorbar(im,cax=fig.add_axes([.90,.17,.02,.27]),ticks=[0,25,50,75,100]);cb.ax.tick_params(labelsize=8)
cb.set_label('Exposed population-days (%)',fontsize=8)
save(fig,'Figure_4_temporal_exposure')

# Six-year comparisons retain the distinction between concentration and population-time.
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,4.6),sharey=True)
fig.subplots_adjust(left=.185,right=.965,bottom=.20,top=.85,wspace=.25)
q=C.set_index('area').loc[order];yy=np.arange(8)
a.hlines(yy,np.minimum(q.ghap_six_year_mean,q.acag_six_year_mean),np.maximum(q.ghap_six_year_mean,q.acag_six_year_mean),color='#c1cdd2',lw=2)
a.scatter(q.ghap_six_year_mean,yy,s=30,color=teal,zorder=3)
a.scatter(q.acag_six_year_mean,yy,s=28,marker='s',color=orange,zorder=3)
a.set_yticks(yy,order,fontsize=9);a.set_ylim(7.6,-.6)
a.set_xlabel('Six-year mean PM₂.₅\n(µg/m³)',fontsize=9)
b.barh(yy,q.population_days_above15_annual_weights/1e6,color=teal,height=.53)
for i,v in enumerate(q.population_days_above15_annual_weights/1e6):b.text(v+3,i,f'{v:.1f}',va='center',fontsize=8)
b.set_xlim(0,255);b.set_xlabel('Population-days >15 µg/m³\n(millions)',fontsize=9)
panel(a,'a');panel(b,'b')
fig.legend(handles=product_handles,loc='upper center',bbox_to_anchor=(.54,.99),ncol=2,frameon=False,fontsize=9)
save(fig,'Figure_5_district_exposure')

# Trends, intermediate classifications and spatial ranking tell different stories.
fig=plt.figure(figsize=(7.2,6.5))
a=fig.add_axes([.11,.62,.36,.27]);b=fig.add_axes([.61,.62,.35,.27]);c=fig.add_axes([.20,.12,.73,.34])
q=comp[comp.area=='Midlands Province'].sort_values('year');years=q.year.to_numpy()
for val,color,marker in [(q.ghap_v1_pw_pm25,teal,'o'),(q.comparator_pw_pm25,orange,'s')]:a.plot(years,val,color=color,marker=marker,ms=4,lw=1.1)
a.axhline(5,color=muted,ls='--',lw=.8);a.set_ylim(4,16);a.set_ylabel('Annual PM₂.₅ (µg/m³)',fontsize=9)
gh=annual[annual.area=='Midlands Province'].sort_values('year')
b.plot(years,100*gh.population_above_annual_10.to_numpy()/gh.population.to_numpy(),'o-',color=teal,ms=4,lw=1.1)
b.plot(years,100*q.comparator_population_above_10/q.comparator_population,'s-',color=orange,ms=4,lw=1.1)
b.set_ylim(0,105);b.set_ylabel('Population >10 µg/m³ (%)',fontsize=9)
for ax in [a,b]:ax.set_xticks(years);ax.tick_params(axis='x',labelrotation=45);ax.grid(axis='y',alpha=.15)
s=S.set_index('area').loc[order]
c.hlines(yy,np.minimum(s.ghap_rank,s.acag_rank),np.maximum(s.ghap_rank,s.acag_rank),color='#becbd1',lw=2)
c.scatter(s.ghap_rank,yy-.085,s=32,color=teal,zorder=3);c.scatter(s.acag_rank,yy+.085,s=30,marker='s',color=orange,zorder=3)
c.set_yticks(yy,order,fontsize=9);c.set_ylim(7.6,-.6);c.set_xlim(.6,8.4);c.set_xticks(range(1,9))
c.set_xlabel('Six-year concentration rank (1 = highest)',fontsize=9)
for ax,letter in [(a,'a'),(b,'b'),(c,'c')]:panel(ax,letter)
fig.legend(handles=product_handles,loc='upper center',bbox_to_anchor=(.53,.985),ncol=2,frameon=False,fontsize=9)
save(fig,'Figure_6_product_sensitivity')

# Two large map panels leave room for readable external keys.
fig=plt.figure(figsize=(7.2,5.8));axes=[]
for x,letter in [(.055,'a'),(.535,'b')]:
    ax=fig.add_axes([x,.29,.39,.64]);map_axes(ax,letter,ticks=False);axes.append(ax)
mx=float(pd.read_csv(RES/'cumulative_native_cell_profiles.csv').cumulative_excess_ug_m3_days.max())
im=map_raster(axes[0],RES/'rasters/cumulative_excess_2017_2022.tif','YlOrBr',Normalize(0,np.ceil(mx/1000)*1000))
horizontal_key(fig,im,[.065,.205,.365,.022],'Cumulative excess (µg/m³-days)')
cm=plt.get_cmap('YlGnBu',7)
im=map_raster(axes[1],RES/'rasters/years_above10_2017_2022.tif',cm,BoundaryNorm(np.arange(-.5,7.5),7))
horizontal_key(fig,im,[.545,.205,.365,.022],'Years with annual PM₂.₅ >10 µg/m³',range(7))
north(fig.add_axes([.465,.065,.045,.11]));scale_bar(fig,axes[0],.19,.055,length=50);scale_bar(fig,axes[1],.69,.055,length=50)
save(fig,'Figure_7_cumulative_exposure')

# Intervals replace the main point-only health heatmap.
fig,axs=plt.subplots(1,3,figsize=(7.2,4.9),sharey=True)
fig.subplots_adjust(left=.185,right=.975,bottom=.21,top=.81,wspace=.20)
for ax,letter in zip(axs,'abc'):
    panel(ax,letter);ax.axvline(0,color=muted,ls='--',lw=.7);ax.set_ylim(7.6,-.6)
axs[0].set_yticks(yy,order,fontsize=9)
for source,color,marker,off in [('Agache 2024',teal,'o',-.13),('Lin 2024',orange,'s',.13)]:
    z=H[(H.year==2022)&(H.endpoint=='asthma_ED')&(H.source==source)].set_index('area').loc[order]
    v=100*(z.rr_estimate.to_numpy()-1);lo=100*(z.ci95_lower.to_numpy()-1);hi=100*(z.ci95_upper.to_numpy()-1)
    axs[0].errorbar(v,yy+off,xerr=[v-lo,hi-v],fmt=marker,color=color,ms=3.5,lw=1,capsize=2,label=source)
z=R[(R.period.astype(str)=='2022')&(R.summary=='p95')&(R.function_id=='lin2024_COPD_admissions_w7')].set_index('area').loc[order]
v=100*(z.rr_estimate.to_numpy()-1)
axs[1].errorbar(v,yy,xerr=[v-100*(z.ci95_lower.to_numpy()-1),100*(z.ci95_upper.to_numpy()-1)-v],fmt='o',color=navy,ms=3.5,lw=1,capsize=2)
for product,color,marker,off in [('GHAP',teal,'o',-.13),('ACAG',orange,'s',.13)]:
    z=R[(R.summary=='six_year_mean')&(R.function_id=='orellano2024_IHD_mortality')&(R['product']==product)].set_index('area').loc[order]
    axs[2].hlines(yy+off,100*(z.pi80_lower-1),100*(z.pi80_upper-1),color=color,alpha=.2,lw=4)
    v=100*(z.rr_estimate.to_numpy()-1)
    axs[2].errorbar(v,yy+off,xerr=[v-100*(z.ci95_lower.to_numpy()-1),100*(z.ci95_upper.to_numpy()-1)-v],fmt=marker,color=color,ms=3.5,lw=1,capsize=2,label=product)
for ax,label in zip(axs,['Asthma ED\ncontrast (%)','COPD admissions\ncontrast (%)','IHD mortality\ncontrast (%)']):
    ax.set_xlabel(label,fontsize=9);ax.xaxis.set_major_locator(MaxNLocator(4))
axs[0].legend(loc='lower center',bbox_to_anchor=(.5,1.025),fontsize=7.5,frameon=False,handlelength=1)
axs[2].legend(loc='lower center',bbox_to_anchor=(.5,1.025),fontsize=8,frameon=False,handlelength=1,ncol=1)
save(fig,'Figure_8_health_profiles')

# Separate moving-average fits and single-lag meta-analyses are never joined.
fig,axs=plt.subplots(2,2,figsize=(7.2,6.3));fig.subplots_adjust(left=.115,right=.97,bottom=.16,top=.94,wspace=.34,hspace=.52)
for ax,conditions,letter in [(axs[0,0],['asthma_admissions','asthma_ED'],'a'),(axs[0,1],['COPD_admissions','COPD_ED'],'b')]:
    for off,condition,color in zip([-.08,.08],conditions,[teal,orange]):
        z=F[(F.endpoint==condition)&(F.window_days>0)].sort_values('window_days')
        v=100*(z.rr_per_10.to_numpy()-1)
        ax.errorbar(np.arange(3)+off,v,yerr=[v-100*(z.ci95_lower_per_10.to_numpy()-1),100*(z.ci95_upper_per_10.to_numpy()-1)-v],fmt='o',color=color,ms=4,capsize=3,lw=1)
    panel(ax,letter);ax.grid(axis='x',visible=False);ax.grid(axis='y',alpha=.15)
    ax.set_xticks(range(3),['0','0–2','0–6']);ax.set_xlabel('Averaging window (lag days)',fontsize=8)
    ax.set_ylabel('Excess rate per 10 µg/m³ (%)',fontsize=8);ax.set_ylim(bottom=0)
for ax,endpoint,letter,color in [(axs[1,0],'asthma_ED','c',orange),(axs[1,1],'asthma_admissions','d',teal)]:
    z=M[M.endpoint==endpoint].sort_values('lag_days');v=100*(z.rr.to_numpy()-1)
    ax.errorbar(z.lag_days,v,yerr=[v-100*(z.ci95_lower.to_numpy()-1),100*(z.ci95_upper.to_numpy()-1)-v],fmt='o',color=color,ms=4,capsize=3,lw=1)
    panel(ax,letter);ax.grid(axis='x',visible=False);ax.grid(axis='y',alpha=.15);ax.axhline(0,color=muted,ls='--',lw=.8)
    ax.set_xticks(range(5));ax.set_xlabel('Single lag (days)',fontsize=8);ax.set_ylabel('Relative excess per 10 µg/m³ (%)',fontsize=8)
fig.legend(handles=[Line2D([],[],color=teal,marker='o',ls='',label='Admissions'),Line2D([],[],color=orange,marker='o',ls='',label='ED visits')],
           loc='lower center',bbox_to_anchor=(.55,.01),ncol=2,frameon=False,fontsize=9)
save(fig,'Figure_S3_published_windows')

# Full descriptive rank matrix, with all input choices rather than selected maxima.
ar=pd.read_csv(RES/'annual_district_rank_stability.csv');wr=pd.read_csv(RES/'acute_window_rank_stability.csv')
fig=plt.figure(figsize=(7.2,6.6));cm=plt.get_cmap('YlGnBu_r',8);norm=BoundaryNorm(np.arange(.5,9.5),8)
for product,pos,letter in [('GHAP',[.19,.57,.34,.34],'a'),('ACAG',[.61,.57,.34,.34],'b')]:
    ax=fig.add_axes(pos);panel(ax,letter);ax.grid(False)
    mat=ar[ar['product']==product].pivot(index='area',columns='year',values='rank').loc[order]
    im=ax.imshow(mat,aspect='auto',cmap=cm,norm=norm)
    ax.set_xticks(range(6),mat.columns,rotation=45,fontsize=8);ax.set_yticks(yy,order if letter=='a' else ['']*8,fontsize=8);ax.tick_params(length=0)
    for (i,j),v in np.ndenumerate(mat.to_numpy()):ax.text(j,i,f'{v:.0f}',ha='center',va='center',fontsize=8,color=contrast_text(im,v))
ax=fig.add_axes([.19,.18,.76,.26]);panel(ax,'c');ax.grid(False)
mat=wr.pivot(index='area',columns=['window_days','year'],values='rank').loc[order].sort_index(axis=1)
im=ax.imshow(mat,aspect='auto',cmap=cm,norm=norm)
ax.set_xticks(range(18),[str(y)[2:] for _,y in mat.columns],fontsize=8);ax.set_yticks(yy,order,fontsize=8);ax.tick_params(length=0)
for sep in [5.5,11.5]:ax.axvline(sep,color='white',lw=2)
for center,label in [(2.5,'1 day'),(8.5,'3 days'),(14.5,'7 days')]:ax.text(center,9.0,label,ha='center',va='top',fontsize=8,clip_on=False)
for (i,j),v in np.ndenumerate(mat.to_numpy()):ax.text(j,i,f'{v:.0f}',ha='center',va='center',fontsize=7.5,color=contrast_text(im,v))
horizontal_key(fig,im,[.29,.070,.50,.018],'District rank (1 = greatest exposure contrast)',range(1,9))
save(fig,'Figure_S4_rank_stability')

(ROOT/'qa/figure_design_record.json').write_text(json.dumps(dict(main_figures=8,supplementary_figures=4,
    main_health_display='Point estimates with 95% coefficient intervals; long-term 80% prediction intervals also shown',
    panel_titles='Lowercase letters only',embedded_captions=False,map_furniture='Outside map frames',
    dpi=450,formats=['PNG','TIFF','PDF'],interpolation='None',
    additions=['Matched evidence sensitivity','Six-year product rank comparison','All-year and all-window rank matrix','Monthly exposed population-time']),indent=2),encoding='utf-8')
