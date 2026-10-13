"""Scientific figures with all keys, north arrows and scale bars outside maps."""
from pathlib import Path
import os,json,sys
os.environ.pop('PROJ_LIB',None);os.environ.pop('PROJ_DATA',None)
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import pandas as pd
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap,BoundaryNorm,LogNorm,Normalize,TwoSlopeNorm
from matplotlib.patches import PathPatch,Patch
from matplotlib.path import Path as MPath
from matplotlib.ticker import FuncFormatter,MaxNLocator
from shapely.geometry import shape
from shapely.ops import unary_union
from pyproj import Geod

ROOT=Path(__file__).resolve().parents[1];RES=ROOT/'results';PUB=ROOT/'data/public';OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,
    'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#54636d','text.color':'#24333d',
    'axes.labelcolor':'#24333d','xtick.color':'#54636d','ytick.color':'#54636d','savefig.facecolor':'white',
    'pdf.fonttype':42,'ps.fonttype':42})
navy='#173f55';teal='#16838b';orange='#c7772e';muted='#687b87'
band_colors=['#e6f2f1','#a5d4cc','#f4d991','#e5a66d','#cf6e63','#934c6f']
annual_labels=['≤5','>5–10','>10–15','>15–25','>25–35','>35']
daily_labels=['≤15','>15–25','>25–37.5','>37.5–50','>50–75','>75']
features=json.loads((PUB/'study_administrative_units.geojson').read_text())['features']
bounds={f['properties']['shapeName']:shape(f['geometry']) for f in features}
region=unary_union(list(bounds.values()));x0,y0,x1,y1=region.bounds
extent=(x0-.025,x1+.025,y0-.025,y1+.025)
latmid=(y0+y1)/2
annual=pd.read_csv(RES/'annual_exposure_summary.csv')
cells=pd.read_csv(RES/'native_cell_annual_metrics.csv')
risks=pd.read_csv(RES/'transferred_rr_scenarios.csv')
daily=pd.read_csv(RES/'daily_area_exposure_and_scenarios.csv',parse_dates=['date'])
monthly=pd.read_csv(RES/'monthly_exposure.csv')
bands=pd.read_csv(RES/'annual_who_population_bands.csv')
comp=pd.read_csv(RES/'exposure_product_sensitivity.csv')

def parts(g):return list(g.geoms) if g.geom_type=='MultiPolygon' else [g]
def region_path(g):
    verts=[];codes=[]
    for polygon in parts(g):
        for ring in [polygon.exterior,*polygon.interiors]:
            v=np.asarray(ring.coords);verts.extend(v);codes.extend([MPath.MOVETO]+[MPath.LINETO]*(len(v)-2)+[MPath.CLOSEPOLY])
    return MPath(np.asarray(verts),np.asarray(codes))
def outline(ax,g,color=navy,lw=.7,fill=None):
    patch=PathPatch(region_path(g),facecolor=fill or 'none',edgecolor=color,lw=lw,zorder=5)
    if fill is None:patch.set_path_effects([pe.Stroke(linewidth=lw+.65,foreground='white',alpha=.7),pe.Normal()])
    ax.add_patch(patch)
    return patch
def map_axes(ax,title,ticks=True):
    ax.set_title(title,loc='left',fontweight='bold',pad=11)
    ax.set_xlim(extent[:2]);ax.set_ylim(extent[2:]);ax.set_aspect(1/np.cos(np.deg2rad(latmid)))
    ax.set_facecolor('#f4f6f6')
    for n,g in bounds.items():outline(ax,g,lw=.5)
    if ticks:
        ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:.1f}°E'))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y,pos:f'{abs(y):.1f}°S'))
        ax.tick_params(labelsize=8,length=2)
    else:ax.set_xticks([]);ax.set_yticks([])
    return ax
def read(path):
    with rasterio.open(path) as ds:
        a=ds.read(1).astype(float);a[(a==-9999)|~np.isfinite(a)]=np.nan
        if ds.nodata is not None:a[a==ds.nodata]=np.nan
        b=ds.bounds
        return a,(b.left,b.right,b.bottom,b.top),ds.transform
def map_raster(ax,path,cmap,norm):
    a,e,_=read(path)
    im=ax.imshow(a,extent=e,origin='upper',interpolation='nearest',cmap=cmap,norm=norm,zorder=1)
    im.set_clip_path(PathPatch(region_path(region),transform=ax.transData))
    return im
def north(ax):
    ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)
    ax.annotate('',xy=(.5,.75),xytext=(.5,.12),arrowprops={'arrowstyle':'-|>','color':navy,'lw':1.2})
    ax.text(.5,.94,'N',ha='center',va='top',fontweight='bold',fontsize=9)
def scale_bar(fig,map_ax,x,y,length=20,width=.16):
    # A standalone key axis uses the same horizontal graphic scale as its map.
    dx=Geod(ellps='WGS84').fwd((x0+x1)/2,latmid,90,length*1000)[0]-(x0+x1)/2
    physical_width=map_ax.get_position().width*dx/(extent[1]-extent[0])
    bar=fig.add_axes([x,y,physical_width,.06]);bar.set_xlim(0,length);bar.set_ylim(0,1);bar.axis('off')
    bar.plot([0,length],[.5,.5],color=navy,lw=2)
    for value in [0,length/2,length]:bar.plot([value,value],[.4,.6],color=navy,lw=.8)
    bar.text(0,.1,'0',ha='center',va='top',fontsize=8)
    bar.text(length,.1,f'{length} km',ha='center',va='top',fontsize=8)
def footer(fig,text):
    # Source and interpretation notes belong in the manuscript captions.
    pass
def save(fig,name):
    if '--only' in sys.argv and name!=sys.argv[sys.argv.index('--only')+1]:
        plt.close(fig)
        return
    for ext in ['png','pdf','tif']:
        kwargs={'dpi':450} if ext!='pdf' else {}
        if ext=='tif':kwargs['pil_kwargs']={'compression':'tiff_lzw'}
        # Write to a fresh path then replace; OneDrive can reject w+b on hydrated files.
        target=OUT/f'{name}.{ext}'
        temporary=OUT/f'{name}.new.{ext}'
        fig.savefig(temporary,**kwargs)
        temporary.replace(target)
    plt.close(fig);print('Created',name,flush=True)

