# Portable public-path adapter; original_ file is preserved verbatim.
#!/usr/bin/env python3
"""Teacher report figure supplement: exact archived data, no new simulation.

Visual conventions follow the inspected Greffet (2002), Raman (2014), and
Peurifoy (2018) original figures: geometry before response, shared encodings,
original samples, and method-specific evidence boundaries.
"""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon, Rectangle, FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'reproduced'
OUT.mkdir(parents=True, exist_ok=True)
FIG = OUT / 'teacher_figures'
FIG.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT/'data/history/meep_curated_20261030.json'
M = json.loads(SOURCE.read_text())
FONT = '/System/Library/Fonts/Supplemental/Arial Unicode.ttf'
if Path(FONT).exists():
    font_manager.fontManager.addfont(FONT)
plt.rcParams.update({
    'font.family': (font_manager.FontProperties(fname=FONT).get_name() if Path(FONT).exists() else 'DejaVu Sans'),
    'font.size': 11, 'axes.labelsize': 11, 'axes.titlesize': 12,
    'axes.unicode_minus': False, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.linewidth': .7,
    'axes.edgecolor': '#71808A', 'xtick.labelsize': 10,
    'ytick.labelsize': 10, 'svg.fonttype': 'none',
    'savefig.facecolor': 'white', 'mathtext.fontset': 'dejavusans',
})
NAVY='#315C7C'; TEAL='#159B96'; ORANGE='#DE9850'; GRAY='#71808A'
INK='#243746'; LIGHT='#DDE5EB'; PALE='#F3F7F9'
COLORS=[NAVY,TEAL,'#799C67',ORANGE,'#9B7DA9','#B67868']
CATALOG=[]

def save(fig,name,caption,notes):
    for ext in ('png','svg'):
        fig.savefig(FIG/f'{name}.{ext}',dpi=300,bbox_inches='tight',pad_inches=.14)
    plt.close(fig)
    CATALOG.append({'id':name,'png':str((FIG/f'{name}.png').relative_to(ROOT)),
                    'svg':str((FIG/f'{name}.svg').relative_to(ROOT)),
                    'caption_zh':caption,'evidence_boundary_zh':notes})

def grid(ax):
    ax.grid(color=LIGHT,lw=.55,alpha=.65)
    ax.set_axisbelow(True)

def panel(ax,label,title):
    ax.set_title(f'({label})  {title}',loc='left',pad=13,color=INK,fontweight='bold')

def route():
    fig=plt.figure(figsize=(11.8,5.8))
    ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    ax.text(.02,.975,'从探索到可验证的 8–13 μm 定向热发射',va='top',
            fontsize=18,fontweight='bold',color=INK)
    ax.text(.02,.905,'目标：在可加工钛基结构上，同时获得带内发射强度与角域集中度',
            fontsize=12,color=GRAY)
    xs=[.025,.357,.689]; ys=[.565,.125]; w=.285; h=.265
    stages=[
        ('01','Meep 二维局部规律','占空比 F、±30°方向对照\n8 / 10.5 / 13 μm 离散点','已完成探索',NAVY),
        ('02','六结构几何光学预筛','A0–A3 / B / C 几何比较\n8–13 μm 射线角域指标','已完成预筛',NAVY),
        ('03','A3 三维全波核算','COMSOL 材料损耗与功率审查\n高度 / 填充 / 离散角度比较','已完成部分条件',TEAL),
        ('04','加工与材料输入闭环','激光可达形貌、底填充工艺\n真实 TA2 / 膜层光学常数','下一步核实',ORANGE),
        ('05','连续光谱与角域评价','8–13 μm 加密波长、θ / φ\n无偏振平均、锥内功率与占比','尚需补足',ORANGE),
        ('06','小样本代理辅助选点','H / t 交叉补样与整结构留出\n先验证预测，再约束寻优','已启动试点；待扩展',TEAL),
    ]
    # Scientific sequence flows left-to-right, top then bottom.
    for i,(num,title,detail,status,color) in enumerate(stages):
        x=xs[i%3];y=ys[0 if i<3 else 1]
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012,rounding_size=0.014',
                         facecolor='white',edgecolor=LIGHT,lw=1.1))
        ax.add_patch(Rectangle((x,y+h-.008),w,.008,facecolor=color,edgecolor='none'))
        ax.text(x+.018,y+h-.04,num,color=color,fontsize=12,fontweight='bold',va='top')
        ax.text(x+.06,y+h-.04,title,color=INK,fontsize=13,fontweight='bold',va='top')
        ax.text(x+.019,y+h-.104,detail,color=INK,fontsize=11,va='top',linespacing=1.65)
        ax.text(x+.019,y+.025,status,color=color,fontsize=10.5,fontweight='bold')
    for i in [0,1,3,4]:
        yy=ys[0 if i<3 else 1]+h/2
        ax.add_patch(FancyArrowPatch((xs[i%3]+w+.011,yy),(xs[(i+1)%3]-.012,yy),
                        arrowstyle='-|>',mutation_scale=13,color=GRAY,lw=1.2))
    # wrap from full wave to fabrication: never imply a physical result for pending work.
    ax.plot([.83,.83,.17,.17],[.55,.47,.47,.419],color=GRAY,lw=1.1)
    ax.add_patch(FancyArrowPatch((.17,.445),(.17,.408),arrowstyle='-|>',
                                mutation_scale=13,color=GRAY,lw=1.1))
    ax.text(.5,.463,'由数值候选转入实验约束，再回馈计算选点',ha='center',color=GRAY,fontsize=10.5,
            bbox={'facecolor':'white','edgecolor':'none','pad':3})
    ax.text(.025,.048,'证据边界：二维、射线预筛与三维全波分开解读；当前全波离散点尚不等于完整带内角锥性能。',
            color=GRAY,fontsize=10.3)
    save(fig,'research_route',
         '研究路线与当前证据覆盖。蓝色为已完成探索，青色为已有部分证据或已启动试点，橙色为尚待核实或补足环节。',
         '路线图不作为新仿真结果；六结构射线带内角域不能代替真实材料的三维全波角谱。')

def meep_relationships():
    rows=M['stable_2d_direction_pairs']
    fig,axs=plt.subplots(2,2,figsize=(10.5,8.0))
    fig.subplots_adjust(left=.085,right=.98,bottom=.17,top=.88,wspace=.31,hspace=.47)
    fig.suptitle('Meep 二维斜槽：占空比与方向差的局部证据',fontsize=17,fontweight='bold',color=INK,y=.985)
    fig.text(.5,.93,'裸光学钛 · p 偏振 · 10 个离散条件 · 探索性结果，未达到原正式合格判据',
             ha='center',fontsize=11,color=GRAY)
    ax=axs[0,0];panel(ax,'a','真实周期截面与参数')
    geom=next(r for r in M['stable_2d_geometry'] if r['surface_fill']==.65)
    block=Rectangle((-75,-60),150,60,facecolor='#D0DAE2',edgecolor='none')
    ax.add_patch(block)
    v=np.array(geom['vertices_x_depth'],float);v[:,1]*=-1
    for shift in [-50,0,50]:
        vv=v.copy();vv[:,0]+=shift
        p=Polygon(vv,facecolor='white',edgecolor=NAVY,lw=1.6)
        ax.add_patch(p);p.set_clip_path(block)
    ax.plot([-25,-25],[4,-60],ls='--',lw=.8,color=GRAY)
    ax.plot([25,25],[4,-60],ls='--',lw=.8,color=GRAY)
    ax.annotate('',(-16.25,5),(16.25,5),arrowprops={'arrowstyle':'<->','color':NAVY,'lw':.8})
    ax.text(0,8,'开口 32.5 μm；F = 0.65',ha='center',fontsize=9.8,color=NAVY)
    ax.text(-19,-52,'Ti',color=NAVY,fontsize=12)
    ax.text(5,-18,'空气',ha='center',color=GRAY,fontsize=10)
    ax.annotate('',(0,30),(0,14),arrowprops={'arrowstyle':'->','color':GRAY,'lw':1})
    ax.text(-12,28,'表面法线',fontsize=9.5,color=GRAY,ha='right')
    for sx,col,lab in [(-1,NAVY,'+30°'),(1,TEAL,'−30°')]:
        ax.annotate('',(sx*3,2),(sx*15,25),arrowprops={'arrowstyle':'->','color':col,'lw':1.5})
        ax.text(sx*22,20,lab,ha='center',fontsize=10,color=col)
    ax.set(xlim=(-78,78),ylim=(-63,33),xticks=[-50,-25,0,25,50],yticks=[-60,-40,-20,0],
           xlabel='横向坐标 x  (μm)',ylabel='高度 z  (μm)')
    ax.set_aspect('equal',adjustable='box')
    ax.text(.98,.05,'周期 50 μm\n槽轴长 40 μm\n倾角 30°',transform=ax.transAxes,
            ha='right',va='bottom',fontsize=9.5,color=INK,linespacing=1.5,
            bbox={'facecolor':'white','alpha':.8,'edgecolor':'none','pad':2})
    ax=axs[0,1];panel(ax,'b','占空比：存在局部非单调性')
    duty=sorted([r for r in rows if r['lambda_um']==10.5],key=lambda r:r['F'])
    for key,col,lab in [('plus',NAVY,'+30°'),('minus',TEAL,'−30°')]:
        x=[r['F'] for r in duty];y=[100*r[f'A_vol_{key}'] for r in duty]
        ax.plot(x,y,'o-',lw=1.6,ms=6,color=col,label=lab)
        ax.plot(x,[100*r[f'A_flux_{key}'] for r in duty],'o',mfc='white',mec=col,ms=3.8,mew=.8)
    ax.annotate('F = 0.65\n本组最大方向差',(.65,19.6),(.70,23.3),fontsize=9.5,color=INK,
                arrowprops={'arrowstyle':'->','color':GRAY,'lw':.8})
    ax.set(xlabel='表面占空比 F',ylabel='吸收率 α  (%)',xticks=[.50,.65,.80],ylim=(0,28))
    ax.text(.04,.93,'λ = 10.5 μm',transform=ax.transAxes,color=GRAY,fontsize=10)
    grid(ax);ax.legend(frameon=False,loc='lower right',ncol=2,fontsize=10)
    ax=axs[1,0];panel(ax,'c','波长：三个离散点保留同向差异')
    spec=sorted([r for r in rows if r['F']==.65],key=lambda r:r['lambda_um'])
    for key,col,lab in [('plus',NAVY,'+30°'),('minus',TEAL,'−30°')]:
        ax.plot([r['lambda_um'] for r in spec],[100*r[f'A_vol_{key}'] for r in spec],
                'o-',color=col,lw=1.6,ms=6,label=lab)
        ax.plot([r['lambda_um'] for r in spec],[100*r[f'A_flux_{key}'] for r in spec],
                'o',mfc='white',mec=col,ms=3.8,mew=.8)
    ax.set(xlabel='波长 λ  (μm)',ylabel='吸收率 α  (%)',xticks=[8,10.5,13],ylim=(0,28))
    ax.text(.04,.93,'F = 0.65',transform=ax.transAxes,color=GRAY,fontsize=10)
    grid(ax);ax.legend(frameon=False,loc='lower left',ncol=2,fontsize=10)
    ax=axs[1,1];panel(ax,'d','方向差：独立核算相互吻合')
    x=np.array([r['delta_flux_pp'] for r in rows]);y=np.array([r['delta_vol_pp'] for r in rows])
    ax.plot([0,13],[0,13],color=GRAY,lw=.8,ls='--',label='1:1 参考')
    ax.scatter(x,y,s=48,facecolor=ORANGE,edgecolor='white',lw=.8,zorder=3)
    labs=['F .50 / 10.5','F .65 / 8','F .65 / 10.5','F .65 / 13','F .80 / 10.5']
    offsets=[(8,-5),(-84,13),(-84,-19),(-80,8),(7,-16)]
    for xx,yy,lab,off in zip(x,y,labs,offsets):
        ax.annotate(lab,(xx,yy),xytext=off,textcoords='offset points',fontsize=8.8,color=GRAY)
    maxgap=float(np.max(np.abs(x-y)))
    ax.text(.04,.95,f'两种核算的方向差最大偏离\n{maxgap:.4f} 个百分点',transform=ax.transAxes,
            va='top',fontsize=10.2,color=INK,linespacing=1.55)
    ax.set(xlabel='净通量方向差 Δα  (百分点)',ylabel='体吸收方向差 Δα  (百分点)',
           xlim=(0,13),ylim=(0,13),xticks=[0,4,8,12],yticks=[0,4,8,12])
    grid(ax);ax.legend(frameon=False,loc='lower right',fontsize=9)
    fig.text(.085,.057,'实心点：体吸收；空心点：净通量。线段仅连接原始采样点，不是平滑拟合或连续光谱。',
             color=GRAY,fontsize=10)
    fig.text(.085,.021,'Δα = α(+30°) − α(−30°)；两方向对照不等于半球角域、锥内功率或波束半高宽。',
             color=GRAY,fontsize=10)
    save(fig,'meep_key_relationships',
         '二维裸钛斜槽的几何、占空比、离散波长和独立核算对照。F=0.65仅为本组三点中方向差最大者；所有点均为探索性结果，原失败或未合格记录保留。',
         '十个波长/方向条件对应三个占空比设计；非全波段扫描、非三维微腔、非角锥指标。通量与体吸收方向差相近不替代完整收敛验收。')

def cpc(z):
    s,c=np.sin(np.pi/12),np.cos(np.pi/12);k=20*(1+s)
    return (-s*(k+c*z)+np.sqrt(k*k+2*k*c*z))/(c*c)-10

def ray_screen():
    fig=plt.figure(figsize=(11.8,8.4))
    gs=fig.add_gridspec(2,6,height_ratios=[1.05,1.6],hspace=.44,wspace=.72,
                        left=.07,right=.98,bottom=.20,top=.86)
    fig.suptitle('六结构射线预筛：集中比例与可用输出量需要同时比较',fontsize=17,fontweight='bold',color=INK,y=.98)
    fig.text(.5,.925,'8–13 μm · 523 K Planck 加权 · 无偏振平均 · 几何光学模型，不是三维全波结果',
             ha='center',fontsize=11,color=GRAY)
    shapes=['A0','A1','A2','A3','B','C'];top=float(cpc(90))
    for i,name in enumerate(shapes):
        ax=fig.add_subplot(gs[0,i]);col=COLORS[i]
        if name.startswith('A'):
            z=np.linspace(0,90,129);f=10+(top-10)*z/90 if name=='A3' else cpc(z)
            off={'A0':0,'A1':10,'A2':20,'A3':0}[name];cx=-off+off*(f-10)/(top-10)
            ax.fill_betweenx(z,-40,cx-f,color='#E1E7EB');ax.fill_betweenx(z,cx+f,40,color='#E1E7EB')
            ax.plot(cx-f,z,color=NAVY,lw=1.3);ax.plot(cx+f,z,color=NAVY,lw=1.3)
            ax.plot([-off-10,-off+10],[0,0],color=ORANGE,lw=4,solid_capstyle='butt')
            sub='直壁' if name=='A3' else ['曲壁居中','曲壁偏移10','曲壁偏移20'][i]
            ax.text(0,-16,'H 90 / P 80',ha='center',fontsize=8.9,color=GRAY)
        elif name=='B':
            ax.add_patch(Polygon([[-35,0],[35,0],[-5,40],[-25,40]],facecolor='#E1E7EB',edgecolor=NAVY,lw=1.3))
            ax.plot([35,-5],[0,40],color=ORANGE,lw=4,solid_capstyle='butt');sub='斜台柱'
            ax.text(0,-16,'H 40 / P 100',ha='center',fontsize=8.9,color=GRAY)
        else:
            z=np.array([0,40]);left=np.array([-25,-35]);right=np.array([-5,35])
            ax.fill_betweenx(z,-50,left,color='#E1E7EB');ax.fill_betweenx(z,right,50,color='#E1E7EB')
            ax.plot(left,z,color=NAVY,lw=1.3);ax.plot(right,z,color=NAVY,lw=1.3)
            ax.plot([-25,-5],[0,0],color=ORANGE,lw=4,solid_capstyle='butt');sub='斜方腔'
            ax.text(0,-16,'H 40 / P 100',ha='center',fontsize=8.9,color=GRAY)
        ax.text(0,106,f'{name}',ha='center',color=col,fontweight='bold',fontsize=13)
        ax.text(0,92,sub,ha='center',color=INK,fontsize=9.5)
        ax.set(xlim=(-55,55),ylim=(-24,115));ax.set_aspect('equal');ax.axis('off')
    fig.text(.07,.887,'(a)  冻结几何中心截面（尺寸单位：μm；橙线为 M2 假定活性面）',fontsize=11.4,fontweight='bold',color=INK)
    R=M['six_ray_groups']
    for colidx,cone in [(0,15),(3,30)]:
        ax=fig.add_subplot(gs[1,colidx:colidx+3]);panel(ax,'b' if cone==15 else 'c',f'{cone}° 法向角锥的权衡')
        for name,col in zip(shapes,COLORS):
            rr=[next(r for r in R if r['structure']==name and r['surface']==s) for s in ['M1','M2']]
            x=[100*r[f'C{cone}'] for r in rr];y=[100*r[f'Q{cone}'] for r in rr]
            ax.plot(x,y,color=col,lw=.9,alpha=.6)
            ax.scatter(x[0],y[0],s=45,facecolors='white',edgecolors=col,linewidths=1.3,zorder=4)
            ax.scatter(x[1],y[1],s=55,facecolors=col,edgecolors='white',linewidths=.7,zorder=5)
            # Label each endpoint only for M2, plus an explicit M1 cluster note.
            off={'A0':(7,6),'A1':(7,-1),'A2':(7,-6),'A3':(7,0),'B':(7,-6),'C':(7,1)}[name]
            if cone==30:off={'A0':(-26,10),'A1':(8,-2),'A2':(-25,-16),'A3':(7,0),'B':(7,0),'C':(7,0)}[name]
            ax.annotate(name,(x[1],y[1]),xytext=off,textcoords='offset points',fontsize=10,color=col,fontweight='bold')
        ax.axvline(100*np.sin(np.deg2rad(cone))**2,color=GRAY,ls=':',lw=.9)
        ax.set(xlabel=f'角锥内功率占比 C{cone}  (%)',ylabel=f'角锥内绝对归一化功率 Q{cone}  (%)')
        ax.set_xlim((3.5,24) if cone==15 else (17,57));ax.set_ylim((0,4.8) if cone==15 else (0,10.3))
        grid(ax)
        if cone==15:
            ax.text(.43,.92,'A0：15° 内集中更强',transform=ax.transAxes,color=NAVY,fontsize=10.5)
        else:
            ax.text(.40,.92,'A3：30° 内输出更高',transform=ax.transAxes,color=ORANGE,fontsize=10.5)
    handles=[Line2D([],[],marker='o',mfc='white',mec=GRAY,ls='none',label='M1：全表面光学钛'),
             Line2D([],[],marker='o',color=GRAY,ls='none',label='M2：指定活性面 α = 0.90'),
             Line2D([],[],color=GRAY,ls=':',label='C = sin²β：角度无关发射参考')]
    fig.legend(handles=handles,ncol=3,frameon=False,loc='lower center',bbox_to_anchor=(.52,.092),fontsize=10)
    fig.text(.07,.04,'M2 是灰体假设，不代表已实现膜层。A0/A3端口与活性面积匹配；B/C还改变晶格和活性面积，不能只归因为形状。',
             color=GRAY,fontsize=9.8)
    fig.text(.07,.011,'Qβ = 锥内辐射功率 / 同投影面积黑体半球功率；Cβ = 锥内辐射功率 / 本结构半球功率。',
             color=GRAY,fontsize=9.8)
    save(fig,'structure_screen',
         '六结构冻结几何及射线角域权衡。空心点为裸光学钛M1，实心点为假定灰体活性面M2；同一结构两点连线展示表面假设影响。A0偏向15°集中，A3在30°角锥下输出更高。',
         '射线近似不包含波动衍射、干涉、近场、粗糙或真实膜层转移矩阵；M2不是已验证材料方案。横截面不能显示完整三维面积。')

if __name__=='__main__':
    route();meep_relationships();ray_screen()
    metadata={'source':str(SOURCE.relative_to(ROOT)),
              'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              'visual_precedent':'data/report_figures/visual_reference_notes.json',
              'no_new_solve':True,'catalog':CATALOG,
              'stable_2d_selected_pairs':M['stable_2d_direction_pairs'],
              'ray_selected_groups':M['six_ray_groups']}
    (OUT/'figure_catalog.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'figure_count':len(CATALOG),'files':[r['png'] for r in CATALOG]},ensure_ascii=False))
