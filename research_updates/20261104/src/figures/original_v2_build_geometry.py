from pathlib import Path
import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc, FancyArrowPatch
from matplotlib.collections import PatchCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).parent/'figures';OUT.mkdir(exist_ok=True)
plan=json.load(open(ROOT/'meep_screen/lab_comsol_ml_remaining4_mem24r2_20261103/plan.json'))
g=plan['geometry']; pitch=float(g['pitch_um']); w=float(g['opening_flat_to_flat_um']); b=float(g['bottom_flat_to_flat_um'])
H=float(plan['cases'][0]['height_um']);t=float(plan['cases'][0]['film_um'])
font_manager.fontManager.addfont('/System/Library/Fonts/Supplemental/Arial Unicode.ttf')
plt.rcParams.update({'font.family':'Arial Unicode MS','font.size':12,'axes.unicode_minus':False,'svg.fonttype':'none','text.color':'#172A3A','axes.labelcolor':'#172A3A','figure.facecolor':'white','savefig.facecolor':'white'})
TI='#7893A6'; DARK='#244D65'; FILM='#DE9850'; AIR='#EAF5F8'; GREY='#71808A'; TEAL='#159B96'; RED='#D06851'; BLUE='#3573AA'
def hexv(flat,z=0):
    a=np.pi/6+np.arange(6)*np.pi/3
    return np.c_[flat/np.sqrt(3)*np.cos(a),flat/np.sqrt(3)*np.sin(a),np.full(6,z)]
def dim(ax,a,b,text,off=(0,0),color=DARK,fs=12):
    a=np.array(a)+off;b=np.array(b)+off
    ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'<->','lw':1,'color':color,'shrinkA':0,'shrinkB':0})
    mid=(a+b)/2;ax.text(mid[0],mid[1],text,ha='center',va='center',fontsize=fs,color=color,bbox={'fc':'white','ec':'none','pad':1.5})
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=230,bbox_inches='tight',pad_inches=.13)
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight',pad_inches=.13)
    plt.close(fig)

def cavity3d(ax,height=H,film=t,cut=True):
    top=hexv(w,height); bot=hexv(b,0); outer=hexv(pitch,height)
    ft=hexv(b+(w-b)*film/height,film)
    faces=[];cols=[]
    for j in range(6):
        k=(j+1)%6
        if cut and j in (3,4,5):continue
        faces.append([bot[j],bot[k],top[k],top[j]])
        cols.append(['#93AAB9','#6E899E','#A8BAC5'][j%3])
    ax.add_collection3d(Poly3DCollection(faces,facecolors=cols,edgecolors='#597B91',linewidths=.7,alpha=.72))
    # Upper metal land; front sectors remain visible to outline the real opening.
    lands=[[top[j],top[(j+1)%6],outer[(j+1)%6],outer[j]] for j in range(6)]
    ax.add_collection3d(Poly3DCollection(lands,facecolors='#B9C9D1',edgecolors='#6E899E',linewidths=.6,alpha=.85))
    ax.add_collection3d(Poly3DCollection([ft],facecolors=FILM,edgecolors='#AD6F29',linewidths=.9,alpha=1))
    for j in range(6):
        k=(j+1)%6
        ax.plot(*np.array([bot[j],bot[k]]).T,color=DARK,lw=.7,alpha=.7)
    # The omitted front walls are indicated by fine dashed edges.
    for j in (3,4,5):
        ax.plot(*np.array([bot[j],top[j]]).T,color='#9AB0BF',lw=.9,ls='--',alpha=.8)
    ax.set(xlim=(-48,48),ylim=(-48,48),zlim=(-4,height+9))
    ax.set_box_aspect((96,96,height+13));ax.view_init(elev=22,azim=-65)
    ax.set_axis_off()
    return top,bot,ft

fig=plt.figure(figsize=(10.8,5.5))
fig.text(.025,.965,'A  三维剖切示意',fontsize=14.5,weight='bold')
fig.text(.39,.965,'B  周期单元俯视',fontsize=14.5,weight='bold')
fig.text(.70,.965,'C  中央垂直截面',fontsize=14.5,weight='bold')
ax=fig.add_axes([-.015,.10,.40,.80],projection='3d');top,bot,ft=cavity3d(ax)
ax.text3D(-7,-2,53,'空气腔',fontsize=11,color='#37566C')
ax.text3D(-11,-1,5,'SiO₂',fontsize=11,color='#935E1F')
fig.text(.035,.17,'光学钛壁面',fontsize=12,color=DARK)
fig.text(.035,.105,'去除前侧壁仅用于显示内部',fontsize=10.5,color=GREY)
fig.text(.035,.065,'实际结构仍为完整周期腔',fontsize=10.5,color=GREY)
ax=fig.add_axes([.39,.22,.265,.65]);ax.set_aspect('equal');ax.axis('off')
# Nearest neighbours shown faintly; central hexagon is the true Wigner-Seitz cell.
for ang in np.arange(6)*np.pi/3:
    c=pitch*np.array([np.cos(ang),np.sin(ang)])
    ax.add_patch(Polygon(hexv(w)[:,:2]+c,facecolor='#F5F8FA',edgecolor='#DCE5EB',lw=.7,zorder=0))
ax.add_patch(Polygon(hexv(pitch)[:,:2],facecolor='#D6E1E7',edgecolor=DARK,lw=1.2,ls='--'))
ax.add_patch(Polygon(hexv(w)[:,:2],facecolor=AIR,edgecolor=DARK,lw=1.4))
ax.add_patch(Polygon(hexv(b)[:,:2],facecolor=FILM,edgecolor='#A66E30',lw=.9))
# Dimension across opposite vertical flats; same as source geometry.
dim(ax,(-pitch/2,-53),(pitch/2,-53),'P = 80 μm',fs=11)
ax.plot([-40,-40],[-43,-57],color=GREY,lw=.7);ax.plot([40,40],[-43,-57],color=GREY,lw=.7)
dim(ax,(-w/2,7),(w/2,7),'D = 69.55 μm',fs=10.5)
ax.text(0,43,'周期边界',ha='center',fontsize=11,color=GREY)
ax.set(xlim=(-53,53),ylim=(-62,56))
# Cross-section y=0 has horizontal cavity halfwidth equal to its apothem.
ax=fig.add_axes([.70,.17,.28,.72]);ax.set_aspect('equal');ax.axis('off')
ax.add_patch(Polygon([(-40,-6),(40,-6),(40,H),(-40,H)],fc=TI,ec='none',alpha=.36))
ax.add_patch(Polygon([(-b/2,0),(b/2,0),(w/2,H),(-w/2,H)],fc=AIR,ec=DARK,lw=1.3))
ftw=b+(w-b)*t/H
ax.add_patch(Polygon([(-b/2,0),(b/2,0),(ftw/2,t),(-ftw/2,t)],fc=FILM,ec='#AF7331',lw=.7))
ax.text(0,55,'空气',ha='center',fontsize=11,color='#496D83')
ax.text(35,60,'Ti',ha='center',fontsize=11,color=DARK)
ax.plot([-45,-37],[H,H],color=GREY,lw=.8);ax.plot([-45,-9],[0,0],color=GREY,lw=.8)
ax.annotate('',xy=(-43,H),xytext=(-43,0),arrowprops={'arrowstyle':'<->','color':DARK,'lw':1})
ax.text(-47,H/2,'H = 100 μm',rotation=90,ha='center',va='center',fontsize=11,color=DARK)
dim(ax,(-b/2,-11),(b/2,-11),'d = 15 μm',fs=10)
# A leader labels the true very thin filling without visually coating the walls.
ax.annotate('底部填充\nt = 2 μm',xy=(0,1.1),xytext=(23,18),fontsize=11,ha='center',color='#925C20',arrowprops={'arrowstyle':'-','lw':.9,'color':'#AE7434'})
ax.set(xlim=(-53,54),ylim=(-19,109))
fig.text(.39,.11,'虚线为周期单元；橙色仅在腔底',fontsize=10.5,color=GREY)
fig.text(.39,.065,'代表几何 H100/t2；H80、t3按同一规则变化',fontsize=10.5,color=GREY)
save(fig,'geometry_atlas')

# Two ray diagrams use the source-code wavevector convention, not an ambiguous angle sketch.
fig=plt.figure(figsize=(10.8,6.2))
fig.text(.035,.96,'A  法向入射  θ = 0°',fontsize=15,weight='bold')
fig.text(.50,.96,'B  斜向入射  θ = 30°，φ = 30°',fontsize=15,weight='bold')
for left,theta in [(0.015,0),(.48,30)]:
    ax=fig.add_axes([left,.285,.44,.60]);ax.set_aspect('equal');ax.axis('off')
    # Section in the incident plane, footprint depiction of the same frustum.
    ax.add_patch(Polygon([(-40,-H),(40,-H),(40,0),(-40,0)],fc='#CDDCE5',ec='none'))
    ax.add_patch(Polygon([(-w/2,0),(w/2,0),(b/2,-H),(-b/2,-H)],fc=AIR,ec=DARK,lw=1.25))
    ax.add_patch(Polygon([(-b/2,-H),(b/2,-H),((b+(w-b)*t/H)/2,-H+t),(-(b+(w-b)*t/H)/2,-H+t)],fc=FILM,ec='#AD7434',lw=.6))
    ax.plot([-52,52],[0,0],color=DARK,lw=1.1)
    # Outward normal above plane. Incident theta measured to inward -z as in native code.
    nx = 18 if theta == 0 else 0
    ax.annotate('',xy=(nx,62),xytext=(nx,0),arrowprops={'arrowstyle':'-|>','color':GREY,'lw':1.2,'mutation_scale':13})
    ax.text(nx+3,63,'外法线 +z',fontsize=11,color=GREY)
    L=65;tr=np.radians(theta);start=np.array([-L*np.sin(tr),L*np.cos(tr)]);end=np.array([0.,0.])
    ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'-|>','color':BLUE,'lw':3,'mutation_scale':17})
    if theta:
        # Incident light continues down-right; mark angle to inward normal below interface.
        v=np.array([L*np.sin(tr),-L*np.cos(tr)])
        ax.annotate('',xy=v*.64,xytext=(0,0),arrowprops={'arrowstyle':'-|>','color':BLUE,'lw':1.2,'linestyle':'--','mutation_scale':11})
        ax.plot([0,0],[0,-55],color=GREY,lw=1,ls='--')
        ax.add_patch(Arc((0,0),60,60,theta1=-90,theta2=-60,color=RED,lw=1.6))
        ax.text(13,-33,'θ=30°',fontsize=11,color=RED)
        ax.text(-43,40,'入射光 k',color=BLUE,fontsize=12)
        ax.text(-49,-121,'θ 是光线相对法线的倾斜角',fontsize=11,color=DARK)
    else:
        ax.text(-40,40,'入射光 k',color=BLUE,fontsize=12)
        ax.text(-45,-121,'沿着表面法线照入同一腔体',fontsize=11,color=DARK)
    ax.text(0,-57,'同一 A3 腔体',ha='center',fontsize=11,color=DARK)
    ax.set(xlim=(-58,58),ylim=(-126,73))
# Compact azimuth inset explicitly shows horizontal k projection.
ax=fig.add_axes([.795,.06,.155,.21]);ax.set_aspect('equal');ax.axis('off')
ax.add_patch(Polygon(hexv(pitch)[:,:2]/1.5,fc='#EDF4F7',ec='#9EB3C0',lw=.9))
ax.annotate('',xy=(38,0),xytext=(0,0),arrowprops={'arrowstyle':'->','color':GREY,'lw':1})
ax.annotate('',xy=(0,36),xytext=(0,0),arrowprops={'arrowstyle':'->','color':GREY,'lw':1})
ax.text(39,-2,'x',fontsize=10);ax.text(-3,39,'y',fontsize=10)
ang=np.pi/6
ax.annotate('',xy=(38*np.cos(ang),38*np.sin(ang)),xytext=(0,0),arrowprops={'arrowstyle':'-|>','color':BLUE,'lw':2,'mutation_scale':13})
ax.add_patch(Arc((0,0),35,35,theta1=0,theta2=30,color=RED,lw=1.3))
ax.text(17,7,'φ',color=RED,fontsize=12)
ax.set(xlim=(-35,47),ylim=(-38,47))
fig.text(.045,.16,'改变的是照射方向，几何形貌保持相同。',fontsize=12.5,weight='bold')
fig.text(.045,.105,'θ=0° 与 θ=30°是两个离散方向，不是一条完整角分布。',fontsize=11.5,color=GREY)
fig.text(.045,.06,'截面为方向示意；虚线仅定义夹角，不表示实际传播或场解。',fontsize=10,color=GREY)
fig.text(.50,.19,'方位角 φ：俯视时，入射波矢',fontsize=11.5,color=DARK)
fig.text(.50,.145,'水平投影与 +x 的夹角。',fontsize=11.5,color=DARK)
fig.text(.50,.085,'本图按实际求解代码的 k 方向约定绘制。',fontsize=10,color=GREY)
save(fig,'incidence_angles')
meta={'source_plan':str((ROOT/'meep_screen/lab_comsol_ml_remaining4_mem24r2_20261103/plan.json').relative_to(ROOT)), 'geometry':{'pitch_um':pitch,'opening_across_flats_um':w,'bottom_across_flats_um':b,'height_um':H,'floorfill_um':t},'hexagon_vertex_angle_rad':'pi/6+j*pi/3','wavevector_convention':'kin=(k sin(theta) cos(phi), k sin(theta) sin(phi), -k cos(theta)); theta to inward -z, phi of horizontal kin projection','reciprocal_outgoing_direction':'kout=-kin; polar theta same, azimuth phi+180deg','geometric_visualization':'schematic based on frozen source; not COMSOL field rendering; cutaway removes front walls for view only','film':'bottom floorfill only; not conformal coating'}
(OUT/'geometry_evidence.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2))
print('Geometry and incidence figures saved:',OUT)
