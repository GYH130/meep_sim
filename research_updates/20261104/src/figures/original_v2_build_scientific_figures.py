"""Replot completed simulation evidence without manufacturing continuous responses."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from collections import Counter

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
OUT=BASE/'figures';OUT.mkdir(parents=True,exist_ok=True)
font_manager.fontManager.addfont('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
plt.rcParams.update({'font.family':'Arial Unicode MS','font.size':12,
 'axes.titlesize':13,'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,
 'text.color':'#213547','axes.labelcolor':'#213547','xtick.color':'#43515F','ytick.color':'#43515F',
 'axes.edgecolor':'#A5B1BB','axes.linewidth':.7,'axes.spines.top':False,'axes.spines.right':False,
 'axes.titleweight':'bold','axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white',
 'legend.frameon':False,'svg.fonttype':'none'})
NAVY='#285B78';TEAL='#1D968D';ORANGE='#D5843D';GRAY='#8E9EAA';PALE='#EAF0F3';DARK='#213547'
S=json.load(open(ROOT/'outputs/ml_remaining4_mem24r2_20261103/all9_completed_summary_20261104.json'))
R=[r['metrics'] for r in S['individual_results']]
MLDIR=ROOT/'outputs/ml_surrogate_pilot_20261102_update_20261104'
PV=json.load(open(MLDIR/'prospective_prediction_validation.json'))
DS=json.load(open(MLDIR/'dataset.json'))
evidence={'source_summary':'outputs/ml_remaining4_mem24r2_20261103/all9_completed_summary_20261104.json',
 'source_ml':'outputs/ml_surrogate_pilot_20261102_update_20261104',
 'scope':'Exploratory, single 10.5 um wavelength. Straight segments only join measured calculation endpoints; no interpolation claims.',
 'new_solves_started':0}

def save(fig,name):
    for ext in ['png','svg']:
        fig.savefig(OUT/f'{name}.{ext}',dpi=230,bbox_inches='tight',pad_inches=.16)
    plt.close(fig)
def grid(ax,axis='y'):
    ax.grid(axis=axis,color='#E0E7EB',lw=.65);ax.set_axisbelow(True)
def panel(ax,l,title):
    ax.set_title(f'{l}  {title}',loc='left',pad=15)
def get(h,t,theta,pol):
    return next(r for r in R if r['height_um']==h and r['floor_fill_thickness_um']==t and r['theta_deg']==theta and r['polarization']==pol)
def avg(h,theta):return np.mean([get(h,2,theta,p)['A_loss_total'] for p in ['s','p']])*100

# 1. Controlled local comparisons. No uncomputed angle lobes or dense spline curves.
fig,axes=plt.subplots(2,2,figsize=(10.8,8.4))
fig.subplots_adjust(hspace=.64,wspace=.43,bottom=.12,top=.93,left=.09,right=.97)
ax=axes[0,0];panel(ax,'a','固定底填充厚度的高度响应')
hd=[]
for theta,color,marker,label in [(0,NAVY,'o','法向 θ=0°'),(30,TEAL,'s','斜向 θ=30°，φ=30°')]:
    vals=[avg(h,theta) for h in [80,100]]
    ax.plot([80,100],vals,color=color,lw=1.5,marker=marker,ms=8,label=label)
    for h,v in zip([80,100],vals):ax.annotate(f'{v:.2f}',(h,v),xytext=(0,10),textcoords='offset points',ha='center',color=color,fontsize=11)
    ax.text(90,np.mean(vals)-2.1,f'+{vals[1]-vals[0]:.3f} 个百分点',ha='center',color=color,fontsize=10.5)
    hd.append({'theta_deg':theta,'height_um':[80,100],'A_unpolarized_percent':vals,'gain_pp':vals[1]-vals[0]})
ax.set(xlim=(76,104),ylim=(11.5,35),xticks=[80,100],xlabel='腔高 H（μm）',ylabel='非偏振平均吸收率（%）')
ax.legend(loc='upper left',fontsize=10.5);grid(ax)
ax.text(.02,.03,'t = 2 μm',transform=ax.transAxes,color=GRAY,fontsize=10.5)

ax=axes[0,1];panel(ax,'b','固定高度的底填充厚度响应')
th=[]
for key,color,mark,label in [('A_loss_total',NAVY,'o','总吸收'),('A_Ti_H',TEAL,'s','钛表面吸收'),('A_film',ORANGE,'D','SiO₂填充吸收')]:
    vals=[get(80,t,0,'s')[key]*100 for t in [2,3]]
    ax.plot([2,3],vals,color=color,lw=1.5,marker=mark,ms=7,label=label)
    th.append({'component':key,'thickness_um':[2,3],'percent':vals,'gain_pp':vals[1]-vals[0]})
    if key=='A_loss_total':
        for t,v in zip([2,3],vals):ax.annotate(f'{v:.2f}',(t,v),xytext=(0,10),textcoords='offset points',ha='center',color=color,fontsize=11)
ax.set(xlim=(1.85,3.15),ylim=(5,40),xticks=[2,3],xlabel='SiO₂底填充厚度 t（μm）',ylabel='吸收占入射功率（%）')
ax.legend(loc='upper left',fontsize=10.5);grid(ax)
ax.text(.04,.02,'H = 80 μm；法向 s\n增量约94%来自填充层',transform=ax.transAxes,fontsize=10.5,color=GRAY,
        bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':2})

ax=axes[1,0];panel(ax,'c','偏振差异取决于观察方向')
configs=[(80,0),(100,0),(80,30),(100,30)]
ys=[3,2,1,0];dif=[]
for (h,theta),y in zip(configs,ys):
    s,pv=[get(h,2,theta,q)['A_loss_total']*100 for q in ['s','p']]
    ax.plot([s,pv],[y,y],color='#C8D4DB',lw=3,zorder=1)
    ax.scatter(s,y,s=58,color=NAVY,marker='o',zorder=3)
    ax.scatter(pv,y,s=58,facecolors='white',edgecolors=ORANGE,linewidths=1.6,marker='D',zorder=4)
    delta=abs(s-pv);dif.append({'H_um':h,'theta_deg':theta,'s_percent':s,'p_percent':pv,'difference_pp':delta})
    ax.text(32.5,y,f'Δ {delta:.3f}',ha='left',va='center',fontsize=10.5,color=GRAY)
ax.set_yticks(ys,['H80  法向','H100  法向','H80  30°/30°','H100  30°/30°'])
ax.set(xlim=(12.8,38),ylim=(-.6,3.9),xlabel='吸收率（%）；右列Δ单位为百分点')
ax.legend(handles=[Line2D([],[],ls='',marker='o',color=NAVY,label='s'),Line2D([],[],ls='',marker='D',markerfacecolor='white',markeredgecolor=ORANGE,label='p')],loc='upper left',ncols=2,fontsize=10.5)
grid(ax,'x')

ax=axes[1,1];panel(ax,'d','吸收增强与法向对比需要同时评价')
nr=[avg(h,30)/avg(h,0) for h in [80,100]];xx=[avg(h,0) for h in [80,100]]
ax.annotate('',xy=(xx[1],nr[1]),xytext=(xx[0],nr[0]),arrowprops={'arrowstyle':'->','color':GRAY,'lw':1.2,'shrinkA':8,'shrinkB':8})
for h,x,y,c in zip([80,100],xx,nr,[NAVY,TEAL]):
    ax.scatter(x,y,s=110,c=c,zorder=4)
    ax.annotate(f'H{h}/t2\n{x:.2f}%，{y:.3f}',(x,y),xytext=(8,-36) if h==80 else (-8,10),textcoords='offset points',color=c,fontsize=11,ha='left' if h==80 else 'right')
ax.set(xlim=(25.5,31.2),ylim=(.53,.655),xlabel='法向非偏振平均吸收率（%）',ylabel='两方向比值 A₃₀ / A₀')
ax.text(.04,.94,'↑ 相对斜向吸收增加\n→ 法向绝对吸收增加',transform=ax.transAxes,ha='left',va='top',fontsize=10.5,color=GRAY)
grid(ax)
fig.text(.09,.022,'数据均为10.5 μm已完成仿真。细线只连接已算端点，不代表连续拟合或未计算条件。',fontsize=11,color=GRAY)
save(fig,'response_relationships')
evidence['height_response']=hd;evidence['thickness_components']=th;evidence['polarization_pairs']=dif
evidence['two_direction_ratio']=[{'H_um':h,'normal_percent':x,'ratio_A30_A0':y} for h,x,y in zip([80,100],xx,nr)]

# 2. Actual-versus-prior prediction, and honest coverage of independent designs.
fig,(ax,bx)=plt.subplots(1,2,figsize=(10.5,5.6),gridspec_kw={'width_ratios':[1.15,1]})
fig.subplots_adjust(left=.09,right=.98,bottom=.30,top=.89,wspace=.38)
panel(ax,'a','新结构的事先预测与实际计算')
ps=[p for p in PV['predictions'] if p['design_id']=='H100_t2' and p['label']=='A_loss_total']
for method,color,mark,label in [('mean',GRAY,'s','均值基线'),('linear_ridge',TEAL,'^','岭回归'),('gp_matern',ORANGE,'o','GP')]:
    subset=[p for p in ps if p['method']==method]
    actual=np.array([p['observed'] for p in subset])*100;pred=np.array([p['predicted'] for p in subset])*100
    mae=np.mean(abs(actual-pred))
    ax.scatter(actual,pred,color=color,marker=mark,s=65,alpha=.9,label=f'{label}   MAE {mae:.2f} pp',zorder=4)
ax.plot([13.5,35],[13.5,35],color='#AEB9C1',lw=1.2,ls='--',zorder=1)
ax.text(26,24,'预测 = 实际',rotation=38,color=GRAY,fontsize=10.5)
ax.set(xlim=(13.5,35),ylim=(13.5,35),xlabel='实际 A_loss（%）',ylabel='事先预测 A_loss（%）')
ax.set_xticks([15,20,25,30,35]);ax.set_yticks([15,20,25,30,35]);ax.legend(loc='upper left',fontsize=10.5,labelspacing=.65)
grid(ax,'both')
ax.annotate('法向 s/p几乎重合',(29.52,31.06),xytext=(22.5,28.7),arrowprops={'arrowstyle':'-','color':GRAY,'lw':.7},fontsize=10.5,color=GRAY)

panel(bx,'b','独立尺寸样本的覆盖范围')
heights=[60,80,100,120];ts=[1,2,3]
counts=Counter((int(r['height_um']),int(r['film_um'])) for r in DS['baseline_rows'])
for yi,h in enumerate(heights):
    for xi,t in enumerate(ts):
        n=counts[(h,t)]
        if n:
            bx.scatter(t,h,s=380,facecolors=NAVY,edgecolors='white',linewidth=1,zorder=3)
            bx.text(t,h,str(n),ha='center',va='center',color='white',fontsize=12,zorder=4)
        else:bx.scatter(t,h,s=350,facecolors='white',edgecolors='#D6DFE5',lw=1.2,zorder=2)
        if (h,t) in [(80,2),(100,2)]:bx.scatter(t,h,s=615,facecolors='none',edgecolors=ORANGE,lw=2,zorder=5)
bx.set(xlim=(.6,3.4),ylim=(49,131),xticks=ts,yticks=heights,xlabel='底填充厚度 t（μm）',ylabel='腔高 H（μm）')
grid(bx,'both')
fig.text(.57,.105,'实心圆：已有10.5 μm条件；数字为条件数\n橙色环：本批新增的独立尺寸组合',fontsize=10.5,color=GRAY)
fig.text(.09,.025,'仅1个未见过的H100/t2结构、4个方向/偏振条件；不将其当作4个独立结构，也不绘制未经校准的置信带。',fontsize=11,color=GRAY)
save(fig,'ml_evidence')
evidence['prospective_H100_t2']=[{k:p[k] for k in ['case_id','condition','method','predicted','observed','absolute_error_pp']} for p in ps]
evidence['baseline_coverage']=[{'H_um':h,'t_um':t,'n_conditions':n,'new_independent_design_in_batch':(h,t) in [(80,2),(100,2)]} for (h,t),n in sorted(counts.items())]

# 3. Methods at the same condition, plus independent saved-field accounting checks.
fig,axes=plt.subplots(1,3,figsize=(11.3,5.2),gridspec_kw={'width_ratios':[1.8,1,1]})
fig.subplots_adjust(left=.065,right=.98,bottom=.22,top=.85,wspace=.3)
ax=axes[0];panel(ax,'a','同条件的三种功率核算')
y=np.arange(9)[::-1]
for i,(r,yi) in enumerate(zip(R,y)):
    vals=[r[k]*100 for k in ['A_loss_total','A_port','A_flux']]
    ax.plot([min(vals),max(vals)],[yi,yi],color='#CBD6DD',lw=2,zorder=1)
    for v,key,c,marker,fc in zip(vals,['A_loss','A_port','A_flux'],[NAVY,TEAL,ORANGE],['o','D','x'],[NAVY,'none',ORANGE]):
        if marker=='x':ax.scatter(v,yi,s=54,marker=marker,color=c,lw=1.5,zorder=4)
        else:ax.scatter(v,yi,s=50,marker=marker,facecolors=fc,edgecolors=c,lw=1.3,zorder=3)
ax.set_yticks(y,[f'ML{i+1}' for i in range(9)])
ax.set(xlim=(11.5,39.5),ylim=(-.6,8.7),xlabel='吸收占入射功率（%）')
ax.legend(handles=[Line2D([],[],ls='',marker='o',color=NAVY,label='A_loss'),Line2D([],[],ls='',marker='D',markerfacecolor='white',markeredgecolor=TEAL,label='A_port'),Line2D([],[],ls='',marker='x',color=ORANGE,label='原始 A_flux')],loc='upper left',ncols=3,bbox_to_anchor=(0,-.13),fontsize=10,columnspacing=.7,handletextpad=.2)
grid(ax,'x')
for ax,field,title,color,threshold,unit in [
    (axes[1],'max_saved_field_absolute_error','保存场绝对残差',NAVY,.2,'百分点'),
    (axes[2],'max_saved_field_relative_error','保存场相对残差',TEAL,2.,'%')]:
    panel(ax,'b' if field.endswith('absolute_error') else 'c',title)
    vals=np.array([r[field] for r in S['individual_results']])*100
    ax.scatter(vals,y,s=46,color=color,zorder=3)
    ax.axvline(threshold,color=GRAY,lw=1,ls='--')
    ax.set(ylim=(-.6,8.7),xlabel=unit)
    ax.set_yticks(y,[])
    if unit=='百分点':ax.set_xlim(0,.245);ax.set_xticks([0,.1,.2]);tx=.205
    else:ax.set_xlim(0,2.38);ax.set_xticks([0,1,2]);tx=2.04
    ax.text(tx,8.45,f'{threshold:g}{unit}',fontsize=10,color=GRAY,ha='left')
    grid(ax,'x')
fig.text(.065,.035,'九项保存场功率复核均满足所示核算条件；原始净通量警告保留。这里的残差不是网格收敛误差或仿真总不确定度。',fontsize=11,color=GRAY)
save(fig,'numerical_agreement')
evidence['method_values_and_saved_residuals']=[{'case_id':r['case_id'],'A_loss_percent':r['metrics']['A_loss_total']*100,'A_port_percent':r['metrics']['A_port']*100,'A_flux_percent':r['metrics']['A_flux']*100,'saved_abs_pp':r['max_saved_field_absolute_error']*100,'saved_rel_percent':r['max_saved_field_relative_error']*100} for r in S['individual_results']]
(BASE/'figure_data_checks.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2)+'\n')
print('\n'.join(str(OUT/f'{n}.png') for n in ['response_relationships','ml_evidence','numerical_agreement']))
