# Portable public-path adapter; original_ file is preserved verbatim.
from pathlib import Path
import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc, PathPatch
from matplotlib.path import Path as MPath
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'reproduced/v2_figures';OUT.mkdir(parents=True,exist_ok=True)
plan=json.load(open(ROOT/'data/latest_ml/frozen_scientific_plan.json'))
g=plan['geometry']; W=float(g['opening_flat_to_flat_um']); B=float(g['bottom_flat_to_flat_um']); P=float(g['pitch_um']);H=float(plan['cases'][0]['height_um']);T=float(plan['cases'][0]['film_um'])
if Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf').exists():
    font_manager.fontManager.addfont('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
plt.rcParams.update({'font.family':'Arial Unicode MS','font.size':12,'text.color':'#20394B','axes.unicode_minus':False,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white'})
DARK='#244D65';TI='#CCDCE5';AIR='#EFF7FA';FILM='#DC9850';GREY='#7B8C96';BLUE='#347EAC';RED='#CD6E57';TEAL='#209A92'
def save(f,n):
 f.savefig(OUT/(n+'.png'),dpi=220,bbox_inches='tight',pad_inches=.13)
 f.savefig(OUT/(n+'.svg'),bbox_inches='tight',pad_inches=.13);plt.close(f)
def outline(ax,zoff=0,film=False):
 ax.set_aspect('equal');ax.axis('off')
 ax.add_patch(Polygon([(-P/2,-H-6+zoff),(P/2,-H-6+zoff),(P/2,zoff),(-P/2,zoff)],fc=TI,ec='none'))
 ax.add_patch(Polygon([(-W/2,zoff),(W/2,zoff),(B/2,-H+zoff),(-B/2,-H+zoff)],fc=AIR,ec=DARK,lw=1.15))
 if film:
  wt=B+(W-B)*T/H
  ax.add_patch(Polygon([(-B/2,-H+zoff),(B/2,-H+zoff),(wt/2,-H+T+zoff),(-wt/2,-H+T+zoff)],fc=FILM,ec='#A87536',lw=.7))
 return ax

# A physically paired incoming and outgoing mode, with opposite wavevectors.
f=plt.figure(figsize=(10.8,5.4))
f.text(.03,.95,'A  要研究的热发射观察方向',fontsize=15,color=RED)
f.text(.51,.95,'B  数值求解对应的反向入射',fontsize=15,color=BLUE)
for left,outgoing in [(.02,True),(.50,False)]:
 ax=f.add_axes([left,.24,.46,.64]);outline(ax,film=True)
 ax.plot([-50,50],[0,0],color=DARK,lw=1.1)
 ax.annotate('',xy=(0,62),xytext=(0,0),arrowprops={'arrowstyle':'-|>','color':GREY,'lw':1.2,'mutation_scale':12})
 ax.text(3,62,'外法线 +z',fontsize=11,color=GREY)
 angle=np.pi/6;end=(-65*np.sin(angle),65*np.cos(angle))
 a,b=((0,0),end) if outgoing else (end,(0,0))
 ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'-|>','color':RED if outgoing else BLUE,'lw':3.2,'mutation_scale':17})
 ax.add_patch(Arc((0,0),45,45,theta1=90,theta2=120,color=RED,lw=1.5))
 ax.text(-29,25,'30°',fontsize=12,color=RED,ha='center')
 ax.text(-48,66,'发射光  k_out' if outgoing else '试探光  k_in',color=RED if outgoing else BLUE,fontsize=11)
 ax.text(0,-54,'相同几何与材料',ha='center',fontsize=11,color=DARK)
 ax.set(xlim=(-57,57),ylim=(-109,76))
f.text(.04,.20,'实验目标：观察这个方向的辐射有多强',fontsize=12.7,color=RED)
f.text(.52,.20,'数值手段：计算沿反方向照入时吸收多少',fontsize=12.7,color=BLUE)
f.text(.50,.135,'互易材料在适用条件下：这对反向光学通道的发射率 = 吸收率',ha='center',fontsize=12.5,color=DARK)
f.text(.50,.075,'当前θ=0°与30°只是两个取样点；完整角域仍需多个θ、φ方向。',ha='center',fontsize=11.5,color=GREY)
f.text(.50,.02,'代码约定 k_in=(sinθcosφ, sinθsinφ, −cosθ)；本例 φ_in=30° → φ_out=210°，且 k_out=−k_in。',ha='center',fontsize=10.3,color=GREY)
save(f,'reciprocal_angle')

# Fabrication morphology comparison; the orange coatings in B/C are visual examples only.
f=plt.figure(figsize=(10.8,6.0))
f.text(.025,.95,'A  当前计算的理想底填',fontsize=14.1,color=DARK)
f.text(.365,.95,'B  溶液干燥后的可能形貌',fontsize=14.1,color=DARK)
f.text(.705,.95,'C  气相沉积的两类覆膜',fontsize=14.1,color=DARK)
ax=f.add_axes([.025,.24,.30,.66]);outline(ax,film=True)
# Real dimension identifiers, kept to one side.
ax.text(0,7,'开口 D=69.55 μm',ha='center',fontsize=10.5,color=GREY)
ax.annotate('',xy=(-44,0),xytext=(-44,-H),arrowprops={'arrowstyle':'<->','color':GREY,'lw':.85})
ax.text(-49,-H/2,'H=100 μm',rotation=90,ha='center',va='center',fontsize=10.5,color=GREY)
ax.annotate('平整 2 μm\n仅在腔底',xy=(0,-99),xytext=(24,-70),fontsize=11,color='#A86C26',ha='center',arrowprops={'arrowstyle':'-','lw':.8,'color':'#B4803E'})
ax.set(xlim=(-55,55),ylim=(-113,15))
f.text(.03,.20,'不含侧壁膜与顶面膜',fontsize=12,color=DARK)
f.text(.03,.145,'这是本批模型假设，',fontsize=11,color=GREY)
f.text(.03,.10,'不是已经制备出的实物截面。',fontsize=11,color=GREY)
ax=f.add_axes([.365,.24,.30,.66]);outline(ax,film=False)
# Irregular, edge-rising meniscus in the bottom part of the same tapered cavity.
z=np.array([-100,-100,-89,-91,-95,-96,-93,-88]);x=np.array([-7.5,7.5,10.5,7,3,-1,-6,-10.7])
ax.add_patch(Polygon(np.c_[x,z],fc=FILM,ec='#A87536',lw=.8))
# Sparse side-wall residual traces; patch locations follow the true slope.
for zs,ze in [(-78,-65),(-52,-41),(-28,-19)]:
 xx=lambda zz: B/2+(W-B)/2*(zz+H)/H
 ax.plot([-xx(zs),-xx(ze)],[zs,ze],color=FILM,lw=2.4,solid_capstyle='round',alpha=.8)
ax.annotate('弯月面 / 局部堆积',xy=(0,-93),xytext=(0,-58),ha='center',fontsize=10.5,color='#A86C26',arrowprops={'arrowstyle':'-','color':'#B4803E','lw':.8})
ax.annotate('少量侧壁残留',xy=(-22,-46),xytext=(-29,-9),ha='center',fontsize=10,color='#A86C26',arrowprops={'arrowstyle':'-','color':'#B4803E','lw':.8})
ax.set(xlim=(-55,55),ylim=(-113,15))
f.text(.37,.20,'润湿、浓度和干燥过程决定形貌',fontsize=11.4,color=DARK)
f.text(.37,.145,'图中仅示一种可能，',fontsize=11,color=GREY)
f.text(.37,.10,'不代表已预测沉积轮廓。',fontsize=11,color=GREY)
# Two mini-sections preserve the same real cavity aspect ratio; their coating width is exaggerated.
for left,kind in [(.70,'ALD'),(.85,'PVD')]:
 ax=f.add_axes([left,.305,.145,.52]);outline(ax,film=False)
 if kind=='ALD':
  # Nearly conformal film along exposed land, walls and floor, visible width deliberately enlarged.
  ax.plot([-40,-W/2,-B/2,B/2,W/2,40],[0,0,-H,-H,0,0],color=FILM,lw=3.3,solid_joinstyle='round')
 else:
  # Directional flux: strongest at top/upper surface, weaker on deep sidewalls and bottom.
  ax.plot([-40,-W/2],[0,0],color=FILM,lw=4)
  ax.plot([W/2,40],[0,0],color=FILM,lw=4)
  for side in (-1,1):
   for za,zb,lw,alpha in [(0,-25,3.5,.98),(-25,-55,2.0,.8),(-55,-80,1.2,.55),(-80,-100,.65,.35)]:
    x=lambda z: side*(B/2+(W-B)/2*(z+H)/H)
    ax.plot([x(za),x(zb)],[za,zb],color=FILM,lw=lw,alpha=alpha)
  ax.plot([-B/2,B/2],[-H,-H],color=FILM,lw=2.5,alpha=.95)
 ax.text(0,12,kind,ha='center',fontsize=12.5,color=DARK)
 ax.set(xlim=(-45,45),ylim=(-108,17))
f.text(.704,.275,'趋于共形',fontsize=11.4,color='#A86C26')
f.text(.854,.275,'更具方向性',fontsize=11.4,color='#A86C26')
f.text(.705,.20,'PVD定向通量也可达腔底',fontsize=11.4,color=DARK)
f.text(.705,.145,'覆盖随视线与角分布变化；',fontsize=11,color=GREY)
f.text(.705,.10,'均不等于只在腔底填入2 μm。',fontsize=11,color=GREY)
f.text(.50,.028,'所有轮廓采用同一A3尺寸；B、C的薄层及不均匀度为放大概念示意，不能读作实测膜厚或预测形貌。',ha='center',fontsize=10.5,color=GREY)
save(f,'fill_morphologies')
meta={'geometry_source':'data/latest_ml/frozen_scientific_plan.json','geometry_um':{'pitch':P,'H':H,'opening':W,'bottom':B,'ideal_fill':T},'reciprocal_angle':'k_out=-k_in; phi_in30 maps to phi_out210 by actual source-code convention','morphology_scope':'A is the ideal frozen model. B/C are process-concept illustrations, not simulated, measured or predicted deposition profiles. Thin layers exaggerated in B/C for clarity.'}
(OUT/'figure_evidence.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2))
print('saved explanatory figures')
