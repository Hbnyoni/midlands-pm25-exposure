"""Draw the study workflow from the implemented public-data methods.

Four stages read left to right with short labels; colour separates environmental
inputs, processing, external health evidence and outputs. Detail stays in the text.
"""
from pathlib import Path
from io import BytesIO
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Circle
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':9.5,'mathtext.default':'regular',
                     'svg.fonttype':'none','pdf.fonttype':42})
W,H=180,104
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,W),ylim=(0,H));ax.axis('off')
ink='#1F3443'
INPUT=('#E4EEF5','#2F6C8F');PROCESS=('#E7F2EC','#2E7D5B');HEALTH=('#FBEEE2','#C0692B');OUTPUT=('#1F3A4D','#1F3A4D')
PM='PM$_{2.5}$'

def box(x,y,w,h,title,sub=None,style=INPUT,light=False):
    fill,edge=style
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0,rounding_size=2',lw=1,facecolor=fill,edgecolor=edge))
    colour='white' if light else ink
    if sub:
        ax.text(x+w/2,y+h/2+2.4,title,ha='center',va='center',fontweight='bold',fontsize=10,color=colour)
        ax.text(x+w/2,y+h/2-3.0,sub,ha='center',va='center',fontsize=9,color=colour if light else '#3C5262')
    else:ax.text(x+w/2,y+h/2,title,ha='center',va='center',fontweight='bold',fontsize=10,color=colour,linespacing=1.15)
    return x,y,w,h

def arrow(a,b,colour=ink,style='-|>',ls='-',connection='arc3'):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=9,lw=1,color=colour,ls=ls,
                                 connectionstyle=connection,shrinkA=0,shrinkB=0))

def stage(x,w,number,label,colour):
    ax.add_patch(Circle((x+4,98),3.2,facecolor=colour,edgecolor='none'))
    ax.text(x+4,98,str(number),ha='center',va='center',fontsize=9.5,fontweight='bold',color='white')
    ax.text(x+9.5,98,label,ha='left',va='center',fontsize=10.5,fontweight='bold',color=colour)
    ax.plot([x,x+w],[93.2,93.2],color=colour,lw=1.1,solid_capstyle='butt')

# Column geometry (mm).
C1,C2,C3,C4=(4,40),(52,30),(90,38),(136,40)
stage(*C1,1,'Inputs',INPUT[1]);stage(*C2,2,'Integration',PROCESS[1])
stage(*C3,3,'Exposure metrics',PROCESS[1]);stage(*C4,4,'Outputs',OUTPUT[1])

rows=[70,49,28];h=17
inputs=[box(C1[0],rows[0],C1[1],h,PM,'GHAP daily, ACAG annual'),
        box(C1[0],rows[1],C1[1],h,'Population','WorldPop, 2017 to 2022'),
        box(C1[0],rows[2],C1[1],h,'Districts','8 reporting areas')]
overlay=box(C2[0],rows[2],C2[1],rows[0]+h-rows[2],'',style=PROCESS)
mid=overlay[1]+overlay[3]/2
ax.text(C2[0]+C2[1]/2,mid+7,'Grid\noverlay',ha='center',va='center',fontweight='bold',fontsize=10,color=ink,linespacing=1.15)
ax.text(C2[0]+C2[1]/2,mid-8,'Population\nallocated to\nnative cells',ha='center',va='center',
        fontsize=9,color='#3C5262',linespacing=1.2)
metrics=[box(C3[0],rows[0],C3[1],h,'Annual and\ndaily means',style=PROCESS),
         box(C3[0],rows[1],C3[1],h,'WHO guideline\ncomparisons',style=PROCESS),
         box(C3[0],rows[2],C3[1],h,'Cumulative\nexcess',style=PROCESS)]
profile=box(C4[0],49,C4[1],38,'District\nexposure profiles',style=OUTPUT,light=True)
contrast=box(C4[0],28,C4[1],17,'Conditional\nhealth contrasts',style=OUTPUT,light=True)
evidence=box(C3[0],4,C3[1],17,'Health evidence','Published functions',style=HEALTH)

for x,y,w,hh in inputs:arrow((x+w,y+hh/2),(C2[0],y+hh/2))
for x,y,w,hh in metrics:
    arrow((C2[0]+C2[1],y+hh/2),(x,y+hh/2))
    target=profile if y>=rows[1] else contrast
    arrow((x+w,y+hh/2),(C4[0],min(max(y+hh/2,target[1]+4),target[1]+target[3]-4)))
# External evidence enters only at the health contrast step.
arrow((C3[0]+C3[1],12.5),(C4[0]+7,28),colour=HEALTH[1],connection='angle,angleA=0,angleB=90')
ax.text(C4[0]+27.5,15.5,'R verification\nof all outputs',ha='center',va='center',fontsize=9,color='#3C5262',linespacing=1.15)
ax.add_patch(FancyBboxPatch((C4[0]+15,9.5),C4[1]-15,12,boxstyle='round,pad=0,rounding_size=2',lw=.9,
                            facecolor='none',edgecolor='#7C8E99',ls=(0,(3,2))))
ax.text(C1[0],12.5,'Environmental data',ha='left',va='center',fontsize=9,color=INPUT[1])
ax.add_patch(FancyBboxPatch((C1[0]+30.5,10.5),4,4,boxstyle='round,pad=0,rounding_size=.8',facecolor=INPUT[0],edgecolor=INPUT[1],lw=.8))
ax.text(C1[0],5.5,'External evidence',ha='left',va='center',fontsize=9,color=HEALTH[1])
ax.add_patch(FancyBboxPatch((C1[0]+30.5,3.5),4,4,boxstyle='round,pad=0,rounding_size=.8',facecolor=HEALTH[0],edgecolor=HEALTH[1],lw=.8))

fig.canvas.draw()
for ext in ['png','tif','pdf','svg']:
    path=OUT/f'Figure_1_workflow.{ext}'
    if ext=='tif':
        buffer=BytesIO();fig.savefig(buffer,format='png',dpi=600,facecolor='white')
        buffer.seek(0)
        with Image.open(buffer) as im:rgb=im.convert('RGB')
        rgb.save(path,compression='tiff_lzw',dpi=(600,600))
    else:fig.savefig(path,dpi=600,facecolor='white')
plt.close(fig)
(ROOT/'qa/workflow_design_record.json').write_text(json.dumps({
    'width_mm':W,'height_mm':H,'raster_dpi':600,'source':'Implemented methods and source manifests',
    'main_figure':1,'embedded_title':False,'embedded_caption':False,'layout':'four stage left to right flow',
    'external_health_evidence_is_separate_input':True,
    'no_fitted_local_health_effects_or_new_statistical_results':True},indent=2),encoding='utf-8')
print('Created workflow in PNG, RGB TIFF, PDF and SVG.')
