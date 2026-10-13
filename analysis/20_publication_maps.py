"""Publication cartography at 180 mm with native, unsmoothed raster support."""
from figure_style import *
from matplotlib.ticker import LogFormatterMathtext
import matplotlib.text
from PIL import Image
from pyproj import Transformer
from shapely.ops import transform as reproject
from matplotlib.patches import Polygon as MPolygon,Rectangle

plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':10,'axes.labelsize':10,
    'axes.titlesize':12,'text.color':'#202B33','axes.labelcolor':'#202B33',
    'xtick.color':'#35424B','ytick.color':'#35424B','svg.fonttype':'none',
    'axes.linewidth':.8,'xtick.major.width':.8,'ytick.major.width':.8})
WIDTH=180/25.4
annual_norm=Normalize(np.floor(cells.annual_pm25.min()),np.ceil(cells.annual_pm25.max()))
annual_cmap='viridis'
band_colors=['#f1f5f7','#b9d9ea','#5ca7ca','#2a788e','#24516f','#182f4d']
records=[]

def north(ax):
    # Compass-style arrow: split dark and light halves, drawn natively so no third-party artwork is embedded.
    ax.set_axis_off();ax.set_xlim(-1,1);ax.set_ylim(-1.15,1.35);ax.set_aspect('equal')
    ax.add_patch(MPolygon([[0,1],[-.42,-.8],[0,-.45]],closed=True,facecolor=navy,edgecolor=navy,lw=.8))
    ax.add_patch(MPolygon([[0,1],[.42,-.8],[0,-.45]],closed=True,facecolor='white',edgecolor=navy,lw=.8))
    ax.text(0,1.08,'N',ha='center',va='bottom',fontweight='bold',fontsize=10,color=navy)

def scale_bar(fig,map_ax,x,y,length=100,width=.16):
    # Alternating segments at the true horizontal scale of the paired map.
    dx=Geod(ellps='WGS84').fwd((x0+x1)/2,latmid,90,length*1000)[0]-(x0+x1)/2
    physical_width=map_ax.get_position().width*dx/(extent[1]-extent[0])
    bar=fig.add_axes([x,y,physical_width,.06]);bar.set_xlim(0,length);bar.set_ylim(0,1);bar.axis('off')
    step=length/4
    for k in range(4):bar.add_patch(Rectangle((k*step,.45),step,.2,facecolor=navy if k%2==0 else 'white',edgecolor=navy,lw=.7))
    ticks=[(0,'0'),(length/2,f'{length/2:g}'),(length,f'{length:g} km')]
    if physical_width*fig.get_figwidth()*25.4<24:ticks.pop(1)  # mid label would collide on short bars
    for value,label in ticks:
        bar.text(value,.32,label,ha='center',va='top',fontsize=9,color='#202B33')

def panelmap(fig,pos,letter,ticks=False):
    ax=fig.add_axes(pos);map_axes(ax,letter,ticks=ticks)
    ax.set_facecolor('white')
    for spine in ax.spines.values():spine.set_visible(ticks);spine.set_color('#9DA9AF');spine.set_linewidth(.7)
    ax.set_title(letter,loc='left',fontweight='bold',fontsize=12,pad=10)
    if ticks:ax.tick_params(labelsize=9,length=3,pad=3)
    return ax

def key(fig,im,pos,label,ticks=None):
    cb=fig.colorbar(im,cax=fig.add_axes(pos),orientation='horizontal',ticks=ticks)
    cb.set_label(label,fontsize=10,labelpad=5);cb.ax.tick_params(labelsize=9,length=3,pad=3)
    cb.outline.set_visible(False)
    return cb


DISTRICT_KEY='    '.join(f'{i} {n}' for i,n in enumerate(bounds,1))
# Hand-tuned label anchors keep names inside small districts and clear of neighbours.
NAME_OFFSET={'Gweru':(-.05,.02),'Shurugwi':(.02,-.03),'Chirumhanzu':(.05,.08),'Zvishavane':(0,-.01)}

def district_numbers(ax):
    # Numbers match the key printed beneath multi-panel maps; names would hide the cell pattern.
    for i,name in enumerate(bounds,1):
        q=bounds[name].representative_point()
        ax.text(q.x,q.y,str(i),ha='center',va='center',fontsize=9,color='#172c38',zorder=7,
            bbox=dict(boxstyle='circle,pad=.12',facecolor='white',edgecolor='none',alpha=.92))

def district_names(ax):
    halo=[pe.withStroke(linewidth=2.6,foreground='white',alpha=.95)]
    for name,g in bounds.items():
        q=g.representative_point();dx,dy=NAME_OFFSET.get(name,(0,0))
        ax.text(q.x+dx,q.y+dy,name,ha='center',va='center',fontsize=9.5,color='#172c38',zorder=7,path_effects=halo)

def add_key_band(fig,mm=8):
    # Grow the canvas downward so existing map elements keep their physical size and place.
    fig.canvas.draw();h=fig.get_figheight()*25.4;H=h+mm;to_fig=fig.transFigure.inverted()
    anchors=[(leg,to_fig.transform(leg.get_bbox_to_anchor().get_points())) for leg in fig.legends]
    texts=[(t,t.get_position()) for t in fig.texts]
    fy=lambda y:(y*h+mm)/H
    for ax in fig.axes:
        b=ax.get_position();ax.set_position([b.x0,fy(b.y0),b.width,b.height*h/H])
    fig.set_figheight(H/25.4)
    for leg,pts in anchors:
        (ax0,ay0),(ax1,ay1)=pts;leg.set_bbox_to_anchor((ax0,fy(ay0),ax1-ax0,(ay1-ay0)*h/H),transform=fig.transFigure)
    for t,(x,y) in texts:t.set_position((x,fy(y)))
    fig.text(.5,(mm*.5)/H,DISTRICT_KEY,ha='center',va='center',fontsize=9,color='#35424B')

def publish(fig,name):
    # At 180 mm all text stays at least 9 pt; at 168 mm in Word it exceeds 8 pt.
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_text() and t.get_fontsize()<9:t.set_fontsize(9)
    fig.canvas.draw()
    for ext in ['png','tif','pdf','svg']:
        tmp=OUT/f'{name}.new.{ext}'; target=OUT/f'{name}.{ext}'
        kw={'dpi':600}
        if ext=='tif':kw['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(tmp,**kw)
        if ext=='tif':
            with Image.open(tmp) as im: rgb=im.convert('RGB')
            rgb.save(tmp,compression='tiff_lzw',dpi=(600,600))
        tmp.replace(target)
    records.append(dict(figure=name,width_mm=180,height_mm=fig.get_figheight()*25.4,
                        minimum_font_pt=9,dpi=600,engine='Python Matplotlib native raster cartography'))
    plt.close(fig);print('Created',name,flush=True)

# Locator and population on exposure support.
fig=plt.figure(figsize=(WIDTH,153/25.4))
loc=fig.add_axes([.03,.47,.31,.46]);loc.set_title('a',loc='left',fontsize=12,fontweight='bold',pad=10)
# OpenStreetMap Standard tiles (cached by 24_fetch_osm_locator_tiles.py) in Web Mercator.
to_merc=Transformer.from_crs(4326,3857,always_xy=True).transform
osm=json.loads((PUB/'osm_tiles_z7/osm_tiles_metadata.json').read_text(encoding='utf-8'))
xs=range(osm['x_range'][0],osm['x_range'][1]+1);ys=range(osm['y_range'][0],osm['y_range'][1]+1)
mosaic=np.vstack([np.hstack([np.asarray(Image.open(PUB/f'osm_tiles_z7/7_{x}_{y}.png').convert('RGB')) for x in xs]) for y in ys])
R=20037508.342789244;size=2*R/2**7
loc.imshow(mosaic,extent=(-R+xs[0]*size,-R+(xs[-1]+1)*size,R-(ys[-1]+1)*size,R-ys[0]*size),origin='upper',interpolation='lanczos',zorder=0)
zimbabwe=unary_union([shape(f['geometry']) for f in json.loads((PUB/'zimbabwe_adm1.geojson').read_text())['features']])
outline(loc,reproject(to_merc,zimbabwe),lw=1.1,color='#202B33')
mid=PathPatch(region_path(reproject(to_merc,region)),facecolor='#c0392b',alpha=.55,edgecolor='#7b1e14',lw=1,zorder=6);loc.add_patch(mid)
(lx0,ly0),(lx1,ly1)=to_merc(24.5,-23.5),to_merc(34,-14)
loc.set_xlim(lx0,lx1);loc.set_ylim(ly0,ly1);loc.set_aspect(1)
loc.set_xticks([]);loc.set_yticks([])
for spine in loc.spines.values():spine.set_visible(True);spine.set_color('#9DA9AF');spine.set_linewidth(.7)
loc.text(.5,-.035,'Basemap © OpenStreetMap contributors',transform=loc.transAxes,ha='center',va='top',fontsize=9,color='#55646d')
ax=panelmap(fig,[.43,.26,.49,.68],'b',ticks=True)
count,e,_=read(RES/'rasters/population_on_pm_grid_2022.tif')
with np.load(ROOT/'qa/spatial_support_primary.npz') as z:area=z['Midlands_Province'].reshape(count.shape)/1e6
density=np.divide(count,area,out=np.full_like(count,np.nan),where=area>0)
zero=ax.imshow(np.where(density==0,1,np.nan),extent=e,origin='upper',cmap=ListedColormap(['#e7eaed']),interpolation='nearest')
zero.set_clip_path(PathPatch(region_path(region),transform=ax.transData))
upper=10**np.ceil(np.log10(np.nanmax(density)))
im=ax.imshow(np.ma.masked_where(~np.isfinite(density)|(density<=0),density),extent=e,origin='upper',
             cmap='cividis',norm=LogNorm(.1,upper),interpolation='nearest')
im.set_clip_path(PathPatch(region_path(region),transform=ax.transData))
district_names(ax)
cb=key(fig,im,[.43,.16,.49,.022],'Residents per km²',10.**np.arange(-1,int(np.log10(upper))+1));cb.formatter=LogFormatterMathtext();cb.update_ticks()
fig.legend(handles=[Patch(facecolor='#e7eaed',label='Modelled zero population')],loc='lower left',bbox_to_anchor=(.025,.36),frameon=False,fontsize=9)
north(fig.add_axes([.04,.15,.08,.14]));scale_bar(fig,ax,.16,.17,length=100)
publish(fig,'Figure_1_study_population')

# A shared scale across years; separate units for exceedance frequency.
fig=plt.figure(figsize=(WIDTH,119/25.4));axes=[]
for i,l in enumerate('abc'):axes.append(panelmap(fig,[.025+.33*i,.30,.295,.63],l))
im=map_raster(axes[0],RES/'rasters/annual_pm25_2017.tif',annual_cmap,annual_norm)
map_raster(axes[1],RES/'rasters/annual_pm25_2022.tif',annual_cmap,annual_norm)
maxdays=np.ceil(cells[cells.year==2022].days_above_15.max()/10)*10
ex=map_raster(axes[2],RES/'rasters/days_above_15_2022.tif','YlOrRd',Normalize(0,maxdays))
key(fig,im,[.04,.20,.59,.023],'Annual PM₂.₅ (µg/m³)',[8,10,12,14,16,18,20])
key(fig,ex,[.70,.20,.26,.023],'Days >15 µg/m³',[0,25,50,75,100])
north(fig.add_axes([.46,.025,.055,.11]));scale_bar(fig,axes[0],.20,.027,length=100)
for a in axes:district_numbers(a)
add_key_band(fig)
publish(fig,'Figure_2_exposure_maps')

# Larger category maps with full-width horizontal population comparisons.
fig=plt.figure(figsize=(WIDTH,149/25.4));axes=[]
for x,l in [(.12,'a'),(.58,'b')]:axes.append(panelmap(fig,[x,.43,.31,.51],l))
bn=BoundaryNorm(np.arange(-.5,6.5),6);cm=ListedColormap(band_colors)
map_raster(axes[0],RES/'rasters/annual_who_band_2022.tif',cm,bn)
ac,e,_=read(PUB/'acag_annual_2022.tif');acb=np.where(np.isfinite(ac),np.searchsorted([5,10,15,25,35],ac,side='left'),np.nan)
im=axes[1].imshow(acb,extent=e,origin='upper',interpolation='nearest',cmap=cm,norm=bn)
im.set_clip_path(PathPatch(region_path(region),transform=axes[1].transData))
ba=fig.add_axes([.18,.215,.77,.15]);ba.set_title('c',loc='left',fontsize=12,fontweight='bold',pad=9)
gh=bands[(bands.area=='Midlands Province')&(bands.year==2022)].sort_values('band').population_percent.to_numpy()
with np.load(RES/'population_weights_comparator_2022.npz') as z:aw=z['Midlands_Province']
acshare=np.array([100*aw[acb.ravel()==i].sum()/aw.sum() for i in range(6)]);left=np.zeros(2)
for i,col in enumerate(band_colors):
    v=np.array([gh[i],acshare[i]]);ba.barh([1,0],v,left=left,color=col,height=.54,edgecolor='white',lw=.8)
    for j,val in enumerate(v):
        if val>7:ba.text(left[j]+val/2,1-j,f'{val:.1f}%',ha='center',va='center',fontsize=9.5,color='#182f4d' if i<2 else 'white')
    left+=v
ba.set_yticks([1,0],['GHAP','ACAG']);ba.set_xlim(0,100);ba.set_xticks([0,25,50,75,100]);ba.set_xlabel('Population share (%)',labelpad=4)
ba.tick_params(labelsize=9);ba.spines['left'].set_visible(False);ba.spines['top'].set_visible(False);ba.spines['right'].set_visible(False)
fig.legend(handles=[Patch(facecolor=c,label=l) for c,l in zip(band_colors,annual_labels)],
    loc='lower center',bbox_to_anchor=(.51,.016),ncol=6,title='Annual PM₂.₅ (µg/m³)',
    fontsize=9,title_fontsize=9,frameon=False,handlelength=1.15,columnspacing=1.15)
north(fig.add_axes([.035,.45,.05,.11]));scale_bar(fig,axes[0],.40,.398,length=100)
for a in axes:district_numbers(a)
add_key_band(fig)
publish(fig,'Figure_3_who_product_comparison')

# Cumulative exposure and persistence remain separate environmental quantities.
fig=plt.figure(figsize=(WIDTH,139/25.4));axes=[]
for x,l in [(.04,'a'),(.55,'b')]:axes.append(panelmap(fig,[x,.29,.41,.65],l))
mx=pd.read_csv(RES/'cumulative_native_cell_profiles.csv').cumulative_excess_ug_m3_days.max()
im=map_raster(axes[0],RES/'rasters/cumulative_excess_2017_2022.tif',annual_cmap,Normalize(0,np.ceil(mx/1000)*1000))
key(fig,im,[.05,.205,.39,.022],'Cumulative excess (µg/m³-days)')
im=map_raster(axes[1],RES/'rasters/years_above10_2017_2022.tif',plt.get_cmap('Blues',7),BoundaryNorm(np.arange(-.5,7.5),7))
key(fig,im,[.56,.205,.39,.022],'Years with annual PM₂.₅ >10 µg/m³',range(7))
north(fig.add_axes([.475,.038,.05,.115]));scale_bar(fig,axes[0],.21,.035,length=100);scale_bar(fig,axes[1],.75,.035,length=100)
for a in axes:district_numbers(a)
add_key_band(fig)
publish(fig,'Figure_7_cumulative_exposure')

fig,axs=plt.subplots(2,3,figsize=(WIDTH,183/25.4));fig.subplots_adjust(left=.02,right=.98,bottom=.17,top=.945,wspace=.10,hspace=.15)
for ax,year,l in zip(axs.ravel(),range(2017,2023),'abcdef'):
    map_axes(ax,l,ticks=False);ax.set_facecolor('white')
    for spine in ax.spines.values():spine.set_visible(False)
    ax.set_title(l,loc='left',fontsize=12,fontweight='bold',pad=10)
    im=map_raster(ax,RES/f'rasters/annual_pm25_{year}.tif',annual_cmap,annual_norm)
key(fig,im,[.22,.10,.56,.018],'Annual PM₂.₅ (µg/m³)',[8,10,12,14,16,18,20])
north(fig.add_axes([.065,.025,.05,.10]));scale_bar(fig,axs.ravel()[0],.83,.035,length=100)
for a in axs.ravel():district_numbers(a)
add_key_band(fig)
publish(fig,'Figure_S1_all_annual_maps')
(ROOT/'qa/cartography_design_record.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
