# Portable public-path adapter; original_ file is preserved verbatim.
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse, Arc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'reproduced/v2_figures';OUT.mkdir(parents=True,exist_ok=True)
if Path('/System/Library/Fonts/Supplemental/Arial Unicode.ttf').exists():
    font_manager.fontManager.addfont('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
plt.rcParams.update({'font.family':'Arial Unicode MS','font.size':12.5,'axes.unicode_minus':False,'svg.fonttype':'none','text.color':'#173246','figure.facecolor':'white','savefig.facecolor':'white'})
NAVY='#244D65';GREY='#748591';TEAL='#198D89';ORANGE='#CE8B40'; BLUE='#3B7DB0'
def save(f,n):
 f.savefig(OUT/(n+'.png'),dpi=220,bbox_inches='tight',pad_inches=.12)
 f.savefig(OUT/(n+'.svg'),bbox_inches='tight',pad_inches=.12)
 plt.close(f)
def node(ax,x,y,w,h,title,sub,edge=NAVY):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008,rounding_size=0.009',facecolor='white',edgecolor=edge,lw=1.1))
 ax.text(x+w/2,y+h*.66,title,ha='center',va='center',fontsize=12.5,color=edge)
 ax.text(x+w/2,y+h*.28,sub,ha='center',va='center',fontsize=10.2,color=GREY)
def arr(ax,a,b,c=NAVY):
 ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,lw=1.2,color=c,shrinkA=3,shrinkB=3))
f,ax=plt.subplots(figsize=(10.8,3.95));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
ax.add_patch(FancyBboxPatch((.012,.66),.976,.32,boxstyle='round,pad=0.005,rounding_size=.012',fc='#F0F4F6',ec='none'))
ax.add_patch(FancyBboxPatch((.012,.13),.976,.48,boxstyle='round,pad=0.005,rounding_size=.012',fc='#ECF6F7',ec='none'))
ax.text(.03,.929,'Windows服务器  COMSOL',fontsize=12.5,color=NAVY)
ax.text(.28,.54,'Mac本地  Python脚本与已安装的scikit-learn库',fontsize=12.5,color=NAVY)
node(ax,.045,.72,.205,.15,'已保存的仿真结果','真实求解  10.5 μm')
node(ax,.755,.72,.205,.15,'新候选点独立求解','产生新数据后回填')
arr(ax,(.75,.795),(.256,.795));ax.text(.50,.83,'追加已完成结果',ha='center',fontsize=11,color=GREY)
node(ax,.045,.225,.205,.18,'质量标签与去重','保留失败  警告  审核状态',TEAL)
node(ax,.28,.225,.205,.18,'训练代理模型','均值  岭回归  GP',TEAL)
node(ax,.515,.225,.205,.18,'整套结构留出验证','同一H/t的所有条件共同留出',TEAL)
node(ax,.75,.225,.205,.18,'提出少量候选点','覆盖缺口与预测价值',TEAL)
arr(ax,(.1475,.716),(.1475,.411));
for a,b in [((.253,.315),(.277,.315)),((.488,.315),(.512,.315)),((.723,.315),(.747,.315))]:arr(ax,a,b,TEAL)
arr(ax,(.855,.412),(.855,.713))
ax.text(.50,.052,'当前10.5 μm基线：8个独立H/t结构、26个方向与偏振条件；软件库已经安装，无需另开图形界面。',ha='center',fontsize=11.1,color=GREY)
save(f,'ml_workflow')

# Definition-only visualization. No response lobe or computed angular interpolation.
f=plt.figure(figsize=(10.8,4.9))
f.text(.035,.94,'A  角锥是连续方向区域',fontsize=14.5,weight='bold')
f.text(.56,.94,'B  本批仅计算两个方向',fontsize=14.5,weight='bold')
ax=f.add_axes([.015,.10,.48,.76],projection='3d')
thc=np.radians(30)
th=np.linspace(0,np.pi/2,22);ph=np.linspace(0,2*np.pi,65)
T,P=np.meshgrid(th,ph);xx=np.sin(T)*np.cos(P);yy=np.sin(T)*np.sin(P);zz=np.cos(T)
ax.plot_wireframe(xx,yy,zz,rstride=4,cstride=4,color='#CEDCE4',linewidth=.4,alpha=.75)
T,P=np.meshgrid(np.linspace(0,thc,16),ph)
ax.plot_surface(np.sin(T)*np.cos(P),np.sin(T)*np.sin(P),np.cos(T),color='#4EA8A2',alpha=.65,edgecolor='none',shade=False)
ph=np.linspace(0,2*np.pi,101)
ax.plot(np.sin(thc)*np.cos(ph),np.sin(thc)*np.sin(ph),np.full_like(ph,np.cos(thc)),color=TEAL,lw=1.3)
ax.plot([0,0],[0,0],[0,1.2],color=NAVY,lw=1.1)
ax.quiver(0,0,0,0,0,1.22,color=NAVY,lw=1.1,arrow_length_ratio=.065)
ax.text(0,0,1.25,'法线 +z',fontsize=11,color=NAVY)
ax.plot([0,np.sin(thc)],[0,0],[0,np.cos(thc)],color=TEAL,lw=1.4)
arc=np.linspace(0,thc,45)
ax.plot(.47*np.sin(arc),np.zeros_like(arc),.47*np.cos(arc),color=ORANGE,lw=1.7)
ax.text(.20,-.02,.52,'θc',color='#A96C24',fontsize=13)
f.text(.08,.72,'目标角锥 Ωc',color=TEAL,fontsize=12)
ax.set_box_aspect((1,1,.73));ax.view_init(elev=22,azim=-62);ax.set(xlim=(-1.05,1.05),ylim=(-1.05,1.05),zlim=(0,1.3));ax.set_axis_off()
f.text(.035,.10,'只画积分区域，不表示本结构的实际辐射瓣。',fontsize=11,color=GREY)
# angular domain rectangle makes the coverage gap explicit.
ax=f.add_axes([.59,.25,.37,.57])
ax.set_xlim(-3,93);ax.set_ylim(-10,370)
ax.set_xlabel('极角 θ（°）',fontsize=12);ax.set_ylabel('方位角 φ（°）',fontsize=12)
ax.set_xticks([0,15,30,60,90]);ax.set_yticks([0,90,180,270,360])
ax.spines[['top','right']].set_visible(False)
ax.axvspan(0,30,facecolor='#D8EEEC',alpha=.9)
ax.text(16,300,'角锥积分需覆盖\n0≤θ≤θc\n0≤φ<360°',ha='center',va='center',color=TEAL,fontsize=11.5)
# For incoming convention, normal phi is arbitrary; show an exact theta30phi30 computational sample.
ax.scatter([0,30],[0,30],s=[50,65],color=[NAVY,BLUE],zorder=4)
ax.annotate('法向\nφ无唯一值',xy=(0,0),xytext=(10,100),fontsize=10.5,color=NAVY,arrowprops={'arrowstyle':'-','color':NAVY,'lw':.8})
ax.annotate('θ30° / φ30°',xy=(30,30),xytext=(43,44),fontsize=11,color=BLUE,arrowprops={'arrowstyle':'-','color':BLUE,'lw':.8})
ax.text(.02,-.31,'右图为入射方向标签；出射互易方向按−k映射。\n两个点不能替代整块角域的加权积分。',transform=ax.transAxes,fontsize=10.2,color=GREY)
save(f,'cone_definition')
print('saved ml_workflow and cone_definition')
