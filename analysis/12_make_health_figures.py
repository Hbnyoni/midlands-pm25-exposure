"""Clean scientific profiles; all explanatory titles remain in DOCX captions."""
from figure_style import *
C=pd.read_csv(RES/'cumulative_area_profiles.csv');R=pd.read_csv(RES/'outcome_specific_health_profiles.csv')
F=pd.read_csv(ROOT/'evidence/health_profile_functions.csv')
order=C[C.area!='Midlands Province'].sort_values('cumulative_excess_ug_m3_days',ascending=False).area.tolist()
fig=plt.figure(figsize=(11.5,6.1))
ax1=fig.add_axes([.05,.23,.245,.70]);ax2=fig.add_axes([.355,.23,.245,.70]);ax3=fig.add_axes([.765,.25,.215,.66])
map_axes(ax1,'a');map_axes(ax2,'b')
mx=float(pd.read_csv(RES/'cumulative_native_cell_profiles.csv').cumulative_excess_ug_m3_days.max())
im=map_raster(ax1,RES/'rasters/cumulative_excess_2017_2022.tif','YlOrBr',Normalize(0,np.ceil(mx/1000)*1000))
cb=fig.colorbar(im,cax=fig.add_axes([.067,.16,.21,.023]),orientation='horizontal');cb.set_label('Excess concentration (µg/m³-days)',fontsize=8);cb.ax.tick_params(labelsize=7)
cm=plt.get_cmap('YlGnBu',7);im=map_raster(ax2,RES/'rasters/years_above10_2017_2022.tif',cm,BoundaryNorm(np.arange(-.5,7.5),7))
cb=fig.colorbar(im,cax=fig.add_axes([.37,.16,.21,.023]),orientation='horizontal',ticks=range(7));cb.set_label('Years with annual PM₂.₅ >10 µg/m³',fontsize=8);cb.ax.tick_params(labelsize=7)
north(fig.add_axes([.625,.095,.04,.13]));scale_bar(fig,ax1,.13,.065,length=50);scale_bar(fig,ax2,.43,.065,length=50)
q=C.set_index('area').loc[order];yy=np.arange(len(q))
ax3.barh(yy,q.population_days_above15_annual_weights/1e6,color=teal,height=.62)
ax3.set_yticks(yy,order,fontsize=8);ax3.invert_yaxis();ax3.set_title('c',loc='left',fontweight='bold',pad=11)
ax3.set_xlabel('Population-days >15 µg/m³\n(millions)',fontsize=8);ax3.grid(axis='x',alpha=.18);ax3.set_axisbelow(True)
save(fig,'Figure_7_cumulative_exposure')

acute=[('lin2024_asthma_ED_w7','Asthma\nED'),('lin2024_COPD_admissions_w7','COPD\nadmissions'),('tian2019_pneumonia_admissions_w3','Pneumonia\nadmissions'),('tian2019_cardiovascular_admissions_w1','Cardiovascular\nadmissions')]
chronic=[('orellano2024_IHD_mortality','Ischaemic\nheart disease'),('orellano2024_cerebrovascular_mortality','Cerebro-\nvascular'),('orellano2024_COPD_mortality','COPD'),('orellano2024_ALRI_mortality','ALRI'),('orellano2024_lung_cancer_mortality','Lung\ncancer')]
fig=plt.figure(figsize=(11.5,5.8))
for number,(spec,period,product,stat,pos,cmap) in enumerate([
    (acute,'2022','GHAP','p95',[.145,.24,.31,.66],'YlOrBr'),
    (chronic,'2017-2022','GHAP','six_year_mean',[.54,.24,.43,.66],'YlGnBu')]):
    ax=fig.add_axes(pos);select=R[(R.period.astype(str)==period)&(R['product']==product)&(R.summary==stat)]
    matrix=np.array([[100*(select[(select.area==a)&(select.function_id==ident)].iloc[0].rr_estimate-1) for ident,_ in spec] for a in order])
    im=ax.imshow(matrix,aspect='auto',cmap=cmap,vmin=0,vmax=np.ceil(matrix.max()/2)*2)
    ax.set_xticks(range(len(spec)),[label for _,label in spec],fontsize=8);ax.tick_params(axis='both',length=0,pad=8)
    ax.set_yticks(range(len(order)),order if number==0 else ['']*len(order),fontsize=8)
    ax.set_title('ab'[number],loc='left',fontweight='bold',pad=11)
    for (i,j),val in np.ndenumerate(matrix):
        rgb=np.asarray(im.cmap(im.norm(val))[:3]);linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
        luminance=float(linear@np.array([.2126,.7152,.0722]))
        ax.text(j,i,f'{val:.2f}',ha='center',va='center',fontsize=8,color='white' if luminance<.22 else '#102a38')
    for sp in ax.spines.values():sp.set_visible(False)
    cb=fig.colorbar(im,cax=fig.add_axes([pos[0],.075,pos[2],.022]),orientation='horizontal')
    cb.set_label('100 × (relative contrast − 1), %',fontsize=8);cb.ax.tick_params(labelsize=7)
save(fig,'Figure_8_health_profiles')

# Coefficients at a common increment: these are separate moving-average fits,
# not independent day-specific lag coefficients or a locally fitted lag curve.
fig,axes=plt.subplots(1,2,figsize=(10,4.4));fig.subplots_adjust(left=.1,right=.96,bottom=.21,top=.9,wspace=.32)
for ax,conditions,letter in [(axes[0],['asthma_admissions','asthma_ED'],'a'),(axes[1],['COPD_admissions','COPD_ED'],'b')]:
    for offset,condition,color in zip([-.09,.09],conditions,[teal,orange]):
        q=F[(F.endpoint==condition)&(F.window_days>0)].sort_values('window_days')
        vals=100*(q.rr_per_10.to_numpy()-1);lo=100*(q.ci95_lower_per_10.to_numpy()-1);hi=100*(q.ci95_upper_per_10.to_numpy()-1)
        ax.errorbar(np.arange(3)+offset,vals,yerr=[vals-lo,hi-vals],fmt='o',capsize=3,color=color,label='Admissions' if condition.endswith('admissions') else 'ED visits')
    ax.set_xticks(range(3),['0','0–2','0–6']);ax.set_xlabel('Exposure averaging window (lag days)')
    ax.set_ylabel('Excess rate per 10 µg/m³ (%)');ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.18)
    ax.set_title(letter,loc='left',fontweight='bold',pad=10)
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,-.005))
save(fig,'Figure_S3_published_windows')
