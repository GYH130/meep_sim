# Portable public-path adapter; original_ file is preserved verbatim.
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

BASE=Path(__file__).resolve().parents[2]/'reproduced'
BASE.mkdir(exist_ok=True);ROOT=Path(__file__).resolve().parents[2];OUT=BASE/'v2_figures';OUT.mkdir(parents=True,exist_ok=True)
if Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf').exists():
    font_manager.fontManager.addfont('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
plt.rcParams.update({'font.family':'Arial Unicode MS','font.size':12,'axes.titlesize':13,
 'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,'axes.spines.top':False,
 'axes.spines.right':False,'axes.edgecolor':'#A5B1BB','axes.linewidth':.7,
 'text.color':'#213547','axes.labelcolor':'#213547','xtick.color':'#43515F','ytick.color':'#43515F',
 'axes.unicode_minus':False,'legend.frameon':False,'svg.fonttype':'none'})
NAVY='#285B78';TEAL='#1D968D';ORANGE='#D5843D';GRAY='#8E9EAA'
DATA=json.load(open(ROOT/'data/ml/dataset.json'))
rows=[r for r in DATA['baseline_rows'] if r['polarization']=='s']
fig,axes=plt.subplots(1,3,figsize=(11.3,4.75))
fig.subplots_adjust(left=.065,right=.985,top=.87,bottom=.31,wspace=.36)
evidence={'scope':'10.5um s-polarization baseline rows only; historic and current mesh/method-quality states retained; no physical law or interpolated surface claimed.','new_solves_started':0,'panels':[]}

def draw(ax,points,color,xkey,panel_name,label_offset=(0,9)):
    pts=sorted(points,key=lambda r:r[xkey])
    xs=[r[xkey] for r in pts];ys=[r['A_loss_total']*100 for r in pts]
    ax.plot(xs,ys,color=color,lw=1.1,ls=(0,(3,3)),alpha=.65,zorder=1)
    out=[]
    for r,x,y in zip(pts,xs,ys):
        quality=r['quality_class']
        if r['method_checked_fit_eligible']:
            ax.scatter(x,y,s=60,color=color,edgecolor='white',lw=.7,zorder=4)
        elif r['method_review_failed']:
            ax.scatter(x,y,s=65,color=color,marker='x',lw=1.7,zorder=4)
        else:ax.scatter(x,y,s=64,facecolors='white',edgecolors=color,lw=1.6,zorder=4)
        ax.annotate(f'{y:.2f}',(x,y),xytext=label_offset,textcoords='offset points',ha='center',color=color,fontsize=10.5)
        out.append({'case_id':r['case_id'],'design_id':r['design_id'],'H_um':r['height_um'],'t_um':r['film_um'],'theta_deg':r['theta_deg'],'phi_deg':r['phi_deg'],'A_loss_percent':y,'quality_class':quality,'method_checked':r['method_checked_fit_eligible'],'method_review_failed':r['method_review_failed']})
    evidence['panels'].append({'name':panel_name,'points':out})

for ax,theta,col,label in [(axes[0],0,NAVY,'a  法向的高度响应'),(axes[1],30,TEAL,'b  斜向的高度响应')]:
    ax.set_title(label,loc='left',pad=15)
    for t,c in [(2,col),(3,ORANGE)]:
        pts=[r for r in rows if r['film_um']==t and r['theta_deg']==theta and r['phi_deg']==(0 if theta==0 else 30)]
        draw(ax,pts,c,'height_um',f'height_theta{theta}_t{t}',(0,-19) if theta==30 and t==2 else (0,10))
    ax.set(xlim=(55,125) if theta==0 else (55,105),xticks=[60,80,100,120] if theta==0 else [60,80,100],xlabel='腔高 H（μm）',ylabel='吸收率 A_loss（%）' if theta==0 else '')
    ax.set_ylim((24,40) if theta==0 else (12,22))
    ax.text(.03,.03,'θ = 0°' if theta==0 else 'θ = 30°，φ = 30°',transform=ax.transAxes,color=GRAY,fontsize=10.5)
    ax.grid(axis='y',color='#E0E7EB',lw=.65);ax.set_axisbelow(True)
    ax.legend(handles=[Line2D([],[],color=col,ls='--',lw=1.1,label='t = 2 μm'),Line2D([],[],color=ORANGE,ls='--',lw=1.1,label='t = 3 μm')],fontsize=10.5,loc='upper left',ncols=2,columnspacing=.8,handlelength=1.8)

ax=axes[2];ax.set_title('c  厚度非单调的历史诊断线索',loc='left',pad=15)
pts=sorted([r for r in rows if r['height_um']==120 and r['theta_deg']==0],key=lambda r:r['film_um'])
ax.plot([r['film_um'] for r in pts],[r['A_loss_total']*100 for r in pts],color=GRAY,lw=1.1,ls=(0,(3,3)),alpha=.65,zorder=1)
for r,c in zip(pts,[GRAY,NAVY,ORANGE]):
    # H120/t1,t2,t3 are all unreviewed in the frozen dataset. Do not promote them.
    draw(ax,[r],c,'film_um','thickness_H120_'+str(r['film_um']))
ax.set(xlim=(.7,3.3),ylim=(29,40),xticks=[1,2,3],xlabel='底填充厚度 t（μm）')
ax.grid(axis='y',color='#E0E7EB',lw=.65);ax.set_axisbelow(True)
ax.text(.04,.03,'H = 120 μm；法向 s\n三个点均尚未通过方法审查',transform=ax.transAxes,color=GRAY,fontsize=10.5)
handles=[Line2D([],[],ls='',marker='o',color=NAVY,ms=7,label='功率方法已审查'),Line2D([],[],ls='',marker='o',markerfacecolor='white',markeredgecolor=NAVY,ms=7,label='方法未审查'),Line2D([],[],ls='',marker='x',color=NAVY,ms=7,label='方法审核失败')]
fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.06,.13),ncols=3,fontsize=11,columnspacing=2)
fig.text(.065,.055,'λ = 10.5 μm，s偏振。细虚线仅引导视线；存在历史网格和核算状态差异，空心点与×仅用于诊断对照。',fontsize=10.8,color=GRAY)
fig.text(.065,.012,'H120的厚度非单调尚不能作为物理规律，需要同方法、同网格复核后再判断。',fontsize=10.8,color=GRAY)
for ext in ['png','svg']:fig.savefig(OUT/f'historical_response.{ext}',dpi=230,bbox_inches='tight',pad_inches=.16)
plt.close(fig)
(BASE/'historical_response_data.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n')
print(OUT/'historical_response.png')
