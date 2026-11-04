# Portable public-path adapter; original_ file is preserved verbatim.
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
R=Path(__file__).resolve().parents[2]; O=R/'reproduced'/'teacher_figures'
O.mkdir(parents=True, exist_ok=True)
if Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf').exists():
    font_manager.fontManager.addfont('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
plt.rcParams.update({'font.family':'Arial Unicode MS','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'svg.fonttype':'none','axes.titleweight':'bold','axes.edgecolor':'#9FAEB9','axes.linewidth':.7})
C=['#285B78','#1D968D']; G='#8999A5'
x=json.load(open(R/'data/history/comsol_audit_20261030.json'))
d=json.load(open(R/'data/ml/dataset.json'))['rows']
f,axs=plt.subplots(2,2,figsize=(11.2,8.2));f.subplots_adjust(wspace=.38,hspace=.64,bottom=.13,top=.94,left=.09,right=.97)
for idx,lam in enumerate([10.5,13]):
 a=axs[0,idx];a.set_title(f'{"ab"[idx]}  t3：{lam:g} μm的高度与方向响应',loc='left',pad=14)
 for h,c in zip([60,100],C):
  z=sorted([q for q in x['controlled_height_wavelength_table'] if q['height_um']==h and q['wavelength_um']==lam],key=lambda q:q['theta_deg'])
  a.plot([q['theta_deg'] for q in z],[100*q['A_loss_mean'] for q in z],'--',lw=1,color=c,alpha=.65)
  for q in z:
   good=q['both_saved_field_strict_pass'];a.scatter(q['theta_deg'],100*q['A_loss_mean'],s=65,marker='o' if good else 'x',color=c,linewidth=1.6,zorder=5)
  a.plot([],[],color=c,marker='o',label=f'H{h}')
 a.set(xlabel='倾角 θ（°）；斜向 φ=30°',ylabel='非偏振平均吸收率（%）',xlim=(-5,65),ylim=(0,47),xticks=[0,30,60]);a.grid(axis='y',alpha=.18);a.legend(ncols=2,fontsize=10,loc='upper right')
 a.text(.03,.03,'● 功率方法已审查；× 方法复核未通过',transform=a.transAxes,fontsize=9,color=G)
a=axs[1,0];a.set_title('c  H80/t3：法向p的稀疏波长点',loc='left',pad=14)
z=x['H80_normal_p_spectral_coverage']['points'];la=[q['wavelength_um'] for q in z]
for key,c,mark,label in [('A_loss_total',C[0],'o','材料总吸收'),('A_film','#D5843D','D','SiO₂'),('A_Ti_H',C[1],'s','钛')]:
 a.plot(la,[q[key]*100 for q in z],ls='--',color=c,lw=1,alpha=.6)
 a.scatter(la,[q[key]*100 for q in z],marker=mark,facecolors='white',edgecolors=c,s=55,label=label,zorder=3)
a.set(xlabel='波长 λ（μm）',ylabel='吸收占入射功率（%）',xticks=la,ylim=(0,52));a.grid(axis='y',alpha=.18);a.legend(ncols=3,fontsize=9,loc='upper center');a.text(.03,.03,'空心点：保存场方法待审核；不是带内积分',transform=a.transAxes,fontsize=9,color=G)
a=axs[1,1];a.set_title('d  ML池设计覆盖（86条件）',loc='left',pad=14)
heights=[60,80,100,120];ts=[1,2,3]
for h in heights:
 for t in ts:
  z=[q for q in d if q['height_um']==h and q['film_um']==t]
  if not z:
   a.scatter(h,t,s=90,facecolors='none',edgecolors='#D1DADF',marker='s');continue
  n=len(z);good=sum(q['saved_field_power_method_validated'] is True for q in z)
  a.scatter(h,t,s=450+n*5,c=C[0],alpha=.85)
  a.text(h,t,f'{n}\n({good})',ha='center',va='center',fontsize=9.5,color='white')
a.set(xlabel='腔高 H（μm）',ylabel='底填充厚度 t（μm）',xticks=heights,yticks=ts,xlim=(48,132),ylim=(.55,3.45));a.grid(alpha=.12);a.set_axisbelow(True)
a.text(.02,.02,'圆内：条件数（已审查数）\n空框：未计算尺寸',transform=a.transAxes,fontsize=8.5,color=G,va='bottom',bbox={'facecolor':'white','edgecolor':'none','alpha':.88,'pad':1})
f.text(.09,.027,'所有点为实际完成结果；线段仅连接采样点。方法审核与网格收敛是不同层级，当前尚无正式合格结构。',color=G,fontsize=10)
for ext in ['png','svg']:f.savefig(O/f'historical_comsol_scope.{ext}',dpi=240,bbox_inches='tight',pad_inches=.15)
plt.close(f)
