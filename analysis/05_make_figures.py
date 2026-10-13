"""Clean publication graphics: panel letters only, explanatory text in captions."""
from figure_style import *
TOTAL='Midlands Province';YEAR=2022
from matplotlib.lines import Line2D
prov=annual[annual.area==TOTAL].sort_values('year')
unit_names=[n for n in bounds]
latest=annual[(annual.year==YEAR)&(annual.area!=TOTAL)].set_index('area')
all_values=cells.annual_pm25.to_numpy()
annual_norm=Normalize(max(0,np.floor(all_values.min())),np.ceil(all_values.max()))
def bar_color(fig,im,position,label,ticks=None):
    ca=fig.add_axes(position);cb=fig.colorbar(im,cax=ca,orientation='horizontal')
    cb.set_label(label,fontsize=8);cb.ax.tick_params(labelsize=8)
    if ticks is not None:cb.set_ticks(ticks)
    return cb

# 1. National context and population density on the exposure support.
fig=plt.figure(figsize=(7.2,6.0))
loc=fig.add_axes([.045,.61,.25,.30]);loc.set_title('a',loc='left',fontweight='bold')
countries=json.loads((PUB/'natural_earth_110m_countries.geojson').read_text())['features']
for f in countries:
    g=shape(f['geometry'])
    if not(g.bounds[2]<24 or g.bounds[0]>35 or g.bounds[3]<-24 or g.bounds[1]>-14):
        zim=f['properties'].get('ADMIN')=='Zimbabwe' or f['properties'].get('NAME')=='Zimbabwe'
        outline(loc,g,lw=.5,color='#a6b4bb',fill='#e1eaeb' if zim else '#f5f6f7')
outline(loc,region,lw=.65,color=teal,fill=teal)
loc.set_xlim(24.5,34);loc.set_ylim(-23.5,-14);loc.set_aspect(1);loc.set_xticks([]);loc.set_yticks([])
loc.text(25.8,-16,'Zimbabwe',fontsize=8,zorder=6)
keyax=fig.add_axes([.05,.18,.27,.36]);keyax.axis('off')
for i,name in enumerate(unit_names,1):keyax.text(0,1-(i-1)*.115,f'{i}  {name}',fontsize=8,va='top')
ax=fig.add_axes([.40,.24,.46,.70]);map_axes(ax,'b')
count,e,_=read(RES/f'rasters/population_on_pm_grid_{YEAR}.tif')
with np.load(ROOT/'qa/spatial_support_primary.npz') as z:area=z['Midlands_Province'].reshape(count.shape)/1e6
density=np.divide(count,area,out=np.full_like(count,np.nan),where=area>0)
base=ax.imshow(np.where(density==0,1,np.nan),extent=e,origin='upper',cmap=ListedColormap(['#e8eced']),interpolation='nearest')
base.set_clip_path(PathPatch(region_path(region),transform=ax.transData))
upper=max(1000,10**np.ceil(np.log10(np.nanmax(density))))
im=ax.imshow(np.ma.masked_where(~np.isfinite(density)|(density<=0),density),extent=e,origin='upper',cmap='magma_r',norm=LogNorm(1,upper),interpolation='nearest')
im.set_clip_path(PathPatch(region_path(region),transform=ax.transData))
for i,name in enumerate(unit_names,1):
    q=bounds[name].representative_point()
    ax.text(q.x,q.y,str(i),ha='center',va='center',fontsize=8,color=navy,zorder=6,
            bbox={'boxstyle':'circle,pad=.18','facecolor':'white','edgecolor':'none','alpha':.85})
bar_color(fig,im,[.40,.15,.43,.021],'Residents per km²',ticks=[10**i for i in range(int(np.log10(upper))+1)])
north(fig.add_axes([.90,.38,.05,.14]));scale_bar(fig,ax,.63,.045,length=50)
fig.legend(handles=[Patch(facecolor='#e8eced',label='Modelled zero population')],loc='lower left',bbox_to_anchor=(.045,.11),fontsize=8,frameon=False)
save(fig,'Figure_1_study_population')

# 2. Continuous annual exposure and daily exceedance counts.
fig=plt.figure(figsize=(7.2,5.1));axs=[]
for i,letter in enumerate('abc'):
    ax=fig.add_axes([.045+i*.315,.25,.275,.68]);map_axes(ax,letter,ticks=False);axs.append(ax)
im=map_raster(axs[0],RES/'rasters/annual_pm25_2017.tif','YlOrBr',annual_norm)
map_raster(axs[1],RES/f'rasters/annual_pm25_{YEAR}.tif','YlOrBr',annual_norm)
max_days=int(cells[cells.year==YEAR].days_above_15.max())
exim=map_raster(axs[2],RES/f'rasters/days_above_15_{YEAR}.tif','YlOrRd',Normalize(0,np.ceil(max_days/10)*10))
bar_color(fig,im,[.065,.15,.49,.022],'Annual PM₂.₅ (µg/m³)')
bar_color(fig,exim,[.69,.15,.24,.022],'Days >15 µg/m³')
north(fig.add_axes([.59,.055,.05,.11]));scale_bar(fig,axs[0],.17,.033,length=50)
save(fig,'Figure_2_exposure_maps')

# 3. Exposure categories are compared across source products, with population.
fig=plt.figure(figsize=(7.2,5.0));axs=[]
for x,letter in [(.035,'a'),(.365,'b')]:
    ax=fig.add_axes([x,.28,.285,.65]);map_axes(ax,letter,ticks=False);axs.append(ax)
bn=BoundaryNorm(np.arange(-.5,6.5),6);cm=ListedColormap(band_colors)
map_raster(axs[0],RES/f'rasters/annual_who_band_{YEAR}.tif',cm,bn)
ac,e,_=read(PUB/f'acag_annual_{YEAR}.tif');acb=np.where(np.isfinite(ac),np.searchsorted([5,10,15,25,35],ac,side='left'),np.nan)
im=axs[1].imshow(acb,extent=e,origin='upper',interpolation='nearest',cmap=cm,norm=bn)
im.set_clip_path(PathPatch(region_path(region),transform=axs[1].transData))
ba=fig.add_axes([.75,.36,.20,.49]);ba.set_title('c',loc='left',fontweight='bold')
gh=bands[(bands.area==TOTAL)&(bands.year==YEAR)].sort_values('band').population_percent.to_numpy()
with np.load(RES/f'population_weights_comparator_{YEAR}.npz') as z:acw=z['Midlands_Province']
acshares=np.array([100*acw[acb.ravel()==b].sum()/acw.sum() for b in range(6)])
bottom=np.zeros(2)
for b,color in enumerate(band_colors):
    vals=np.array([gh[b],acshares[b]]);ba.bar([0,1],vals,bottom=bottom,color=color,width=.62,edgecolor='white',lw=.5);bottom+=vals
ba.set_xticks([0,1],['GHAP','ACAG'],fontsize=8);ba.set_ylim(0,100);ba.set_ylabel('Population share (%)',fontsize=8);ba.tick_params(axis='y',labelsize=8)
fig.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(band_colors,annual_labels)],loc='lower center',bbox_to_anchor=(.5,.12),ncol=6,frameon=False,
           title='Annual PM₂.₅ intervals (µg/m³)',fontsize=8,title_fontsize=8,columnspacing=1.2,handlelength=1.3)
north(fig.add_axes([.055,.035,.045,.09]));scale_bar(fig,axs[0],.36,.025,length=50)
save(fig,'Figure_3_who_product_comparison')

# 4. Complete-calendar time series and month-by-year exposure.
fig=plt.figure(figsize=(7.2,5.1));d=daily[daily.area==TOTAL]
ax=fig.add_axes([.09,.57,.84,.34]);ax.set_title('a',loc='left',fontweight='bold')
ax.plot(d.date,d.pw_pm25,lw=.6,color=teal,alpha=.85)
ax.axhline(15,color=orange,lw=1,ls='--',label='15 µg/m³');ax.set_ylabel('PM₂.₅ (µg/m³)')
ax.set_xlim(pd.Timestamp('2017-01-01'),pd.Timestamp('2022-12-31'));ax.legend(frameon=False,fontsize=8,loc='upper left');ax.grid(axis='y',alpha=.12)
hm=fig.add_axes([.09,.15,.67,.28]);hm.set_title('b',loc='left',fontweight='bold')
m=monthly[monthly.area==TOTAL].pivot(index='year',columns='month',values='pm25_mean')
im=hm.imshow(m.to_numpy(),aspect='auto',cmap='YlOrBr',vmin=0,vmax=np.ceil(m.to_numpy().max()/5)*5)
hm.set_xticks(range(12),['J','F','M','A','M','J','J','A','S','O','N','D']);hm.set_yticks(range(6),m.index)
for i in range(6):
    for j in range(12):hm.text(j,i,f'{m.iloc[i,j]:.0f}',ha='center',va='center',fontsize=8,color='white' if m.iloc[i,j]>im.norm.vmax*.65 else navy)
ca=fig.add_axes([.81,.15,.025,.28]);cb=fig.colorbar(im,cax=ca);cb.set_label('PM₂.₅ (µg/m³)',fontsize=8);cb.ax.tick_params(labelsize=8)
save(fig,'Figure_4_temporal_exposure')

# 5. District concentration versus the absolute number of exposed population-days.
order=latest.sort_values('population_days_above_15',ascending=False).index.tolist();pos=np.arange(len(order))
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,4.35),gridspec_kw={'width_ratios':[1,1]},sharey=True)
fig.subplots_adjust(left=.18,right=.97,bottom=.18,top=.86,wspace=.25)
cc=comp[(comp.year==YEAR)&(comp.area!=TOTAL)].set_index('area')
gh=np.array([latest.loc[n,'population_weighted_pm25'] for n in order]);ac=np.array([cc.loc[n,'comparator_pw_pm25'] for n in order])
a.hlines(pos,np.minimum(gh,ac),np.maximum(gh,ac),lw=2,color='#ccd6da');a.scatter(gh,pos,s=21,color=teal);a.scatter(ac,pos,s=21,marker='s',color=orange)
a.set_yticks(pos,order,fontsize=8);a.invert_yaxis();a.set_xlabel('Annual PM₂.₅ (µg/m³)',fontsize=8);a.set_title('a',loc='left',fontweight='bold');a.grid(axis='x',alpha=.12)
b.barh(pos,[latest.loc[n,'population_days_above_15']/1e6 for n in order],height=.58,color=teal,alpha=.85)
b.set_xlabel('Population-days >15 µg/m³\n(millions)',fontsize=8);b.set_title('b',loc='left',fontweight='bold');b.grid(axis='x',alpha=.12)
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=teal,label='GHAP V1'),Line2D([],[],marker='s',ls='',color=orange,label='ACAG V6.GL.02')],
           loc='upper center',bbox_to_anchor=(.58,.98),ncol=2,frameon=False,fontsize=8)
save(fig,'Figure_5_district_exposure')

# 6. The strength of temporal conclusions is tested with another annual product.
fig,axs=plt.subplots(1,2,figsize=(7.2,3.9));fig.subplots_adjust(left=.09,right=.97,bottom=.23,top=.86,wspace=.34)
c=comp[comp.area==TOTAL].sort_values('year');years=c.year.to_numpy()
axs[0].plot(years,c.ghap_v1_pw_pm25,'o-',color=teal,ms=4,label='GHAP V1')
axs[0].plot(years,c.comparator_pw_pm25,'s-',color=orange,ms=4,label='ACAG V6.GL.02')
axs[0].axhline(5,color=muted,ls='--',lw=.8);axs[0].set_ylabel('Annual PM₂.₅ (µg/m³)');axs[0].set_title('a',loc='left',fontweight='bold')
axs[0].legend(frameon=False,fontsize=8);axs[0].set_xticks(years);axs[0].tick_params(axis='x',rotation=45,labelsize=8)
ghpct=100*prov.population_above_annual_10.to_numpy()/prov.population.to_numpy();acpct=100*c.comparator_population_above_10.to_numpy()/c.comparator_population.to_numpy()
axs[1].bar(years-.18,ghpct,width=.34,color=teal);axs[1].bar(years+.18,acpct,width=.34,color=orange)
axs[1].set_ylim(0,100);axs[1].set_ylabel('Population above 10 µg/m³ (%)',fontsize=8);axs[1].set_title('b',loc='left',fontweight='bold');axs[1].set_xticks(years);axs[1].tick_params(axis='x',rotation=45,labelsize=8)
save(fig,'Figure_6_product_sensitivity')

# S1. Full annual spatial sequence, labelled in the manuscript caption.
fig,axs=plt.subplots(2,3,figsize=(7.2,7.3));fig.subplots_adjust(left=.05,right=.95,bottom=.17,top=.96,wspace=.12,hspace=.16)
for ax,year,letter in zip(axs.ravel(),range(2017,2023),'abcdef'):
    map_axes(ax,letter,ticks=False);im=map_raster(ax,RES/f'rasters/annual_pm25_{year}.tif','YlOrBr',annual_norm)
bar_color(fig,im,[.22,.09,.57,.021],'Annual PM₂.₅ (µg/m³)')
north(fig.add_axes([.075,.035,.045,.09]));scale_bar(fig,axs.ravel()[0],.86,.035,length=50)
save(fig,'Figure_S1_all_annual_maps')

# S2. Conditional mortality scenarios are supplementary, with no risk cartogram.
fig=plt.figure(figsize=(7.2,4.8));ax=fig.add_axes([.11,.51,.80,.34])
for ident,color,offset,label in [('orellano2024_respiratory',teal,-.13,'Respiratory mortality'),('orellano2024_allcause',navy,.13,'All-cause mortality')]:
    q=risks[(risks.area==TOTAL)&(risks.window=='annual')&(risks.reference_cap==5)&(risks.function_id==ident)].sort_values('year')
    x=q.year.to_numpy()+offset
    ax.vlines(x,100*(q.pi80_lower-1),100*(q.pi80_upper-1),lw=5,color=color,alpha=.16)
    ax.vlines(x,100*(q.ci95_lower-1),100*(q.ci95_upper-1),lw=1.4,color=color)
    ax.scatter(x,100*(q.rr_estimate-1),s=19,color=color,label=label)
ax.axhline(0,color=muted,lw=.8,ls='--');ax.set_xticks(years);ax.set_ylabel('RR contrast above 1 (%)');ax.set_title('a',loc='left',fontweight='bold')
fig.legend(*ax.get_legend_handles_labels(),loc='upper left',bbox_to_anchor=(.105,.97),ncol=2,frameon=False,fontsize=8)
b=fig.add_axes([.11,.13,.80,.24]);b.set_title('b',loc='left',fontweight='bold')
b.plot(d.date,100*(d.liu2019_respiratory_rr-1),color=teal,lw=.65)
b.fill_between(d.date,100*(d.liu2019_respiratory_ci_lower-1),100*(d.liu2019_respiratory_ci_upper-1),color=teal,alpha=.2,lw=0)
b.set_ylabel('RR contrast above 1 (%)');b.set_xlim(pd.Timestamp('2017-01-01'),pd.Timestamp('2022-12-31'));b.grid(axis='y',alpha=.12)
save(fig,'Figure_S2_supplementary_risk')
(ROOT/'qa/figure_design_record.json').write_text(json.dumps({'main_figures':6,'supplementary_figures':2,
    'panel_titles':'Lowercase panel letters only; single-panel maps have no title',
    'embedded_captions':False,'source_notes':'In manuscript figure captions',
    'map_furniture':'Legends, north arrows and scale bars outside all map frames',
    'grid_display':'Native cells without interpolation; population density aggregated conservatively to exposure support',
    'dpi':450,'formats':['PNG','TIFF','PDF']},indent=2))
