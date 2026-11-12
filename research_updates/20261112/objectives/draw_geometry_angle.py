from pathlib import Path
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon

OUT = Path(__file__).resolve().parent
font_path = Path(os.environ.get('MEEP_CJK_FONT', 'Arial Unicode.ttf'))
if not font_path.is_file():
    raise SystemExit('Set MEEP_CJK_FONT to a CJK font file to regenerate the Chinese diagram.')
font_manager.fontManager.addfont(str(font_path))
font_family = font_manager.FontProperties(fname=str(font_path)).get_name()
plt.rcParams.update({'font.family': font_family, 'font.size': 11,
                     'axes.unicode_minus': False, 'pdf.fonttype': 42})
navy, blue, orange, grey = '#203D59', '#238FA9', '#D98C38', '#D8DFE4'
fig = plt.figure(figsize=(12, 5.6), facecolor='white')
ax = fig.add_axes([.045, .16, .47, .72])
H, b, D, t = 80., 15., 69.5455842631299, 3.
dtop = b + (D-b)*t/H
ax.add_patch(Polygon([(-45,-4), (45,-4), (45,H), (D/2,H),
                     (b/2,0), (-b/2,0), (-D/2,H), (-45,H)],
                    facecolor=grey, edgecolor=navy, linewidth=1.2))
ax.add_patch(Polygon([(-b/2,0), (b/2,0), (dtop/2,t), (-dtop/2,t)],
                    facecolor=blue, edgecolor=blue))
ax.plot([-D/2,-b/2,b/2,D/2], [H,0,0,H], color=navy, lw=1.6)
ax.annotate('', xy=(-D/2,87), xytext=(D/2,87),
            arrowprops={'arrowstyle':'<->','color':navy})
ax.text(0,90, '开口 D ≈ 69.55 μm', ha='center', color=navy)
ax.annotate('', xy=(50,0), xytext=(50,H),
            arrowprops={'arrowstyle':'<->','color':navy})
ax.text(53,40,'腔深 H\n80 μm', va='center', color=navy)
ax.annotate('', xy=(-b/2,-10), xytext=(b/2,-10),
            arrowprops={'arrowstyle':'<->','color':navy})
ax.text(0,-15,'底部对边距 b = 15 μm',ha='center',color=navy)
ax.annotate('SiO₂ 底填充\nt = 3 μm', xy=(-3,1.5), xytext=(-50,25),
            arrowprops={'arrowstyle':'-','color':blue}, color=blue, ha='left')
ax.text(29,35,'Ti',color=navy,fontsize=13)
ax.text(0,54,'空气腔',ha='center',color='#687B89')
ax.set(xlim=(-58,72),ylim=(-22,99),aspect='equal')
ax.axis('off')
fig.text(.065,.93,'(a) 当前六棱台腔体：中央截面',fontsize=14,color=navy)
fig.text(.065,.08,'示例为原基准尺寸；t 是底填充高度，不是侧壁镀膜厚度。',fontsize=10,color='#556877')

ax3 = fig.add_axes([.55,.16,.41,.72],projection='3d')
phi = np.linspace(0,2*np.pi,121)
z = np.linspace(0,1,35)
P,Z = np.meshgrid(phi,z)
for angle,color,alpha in [(30,orange,.15),(15,blue,.3)]:
    R=Z*np.tan(np.deg2rad(angle))
    ax3.plot_surface(R*np.cos(P),R*np.sin(P),Z,color=color,alpha=alpha,
                     linewidth=0,shade=False)
    rr=np.tan(np.deg2rad(angle))
    ax3.plot(rr*np.cos(phi),rr*np.sin(phi),np.ones_like(phi),color=color,lw=1.4)
    ax3.plot([0,rr],[0,0],[0,1],color=color,lw=1.8)
    ax3.text(rr+.025,0,1.02,f'{angle}°',color=color,fontsize=12)
ax3.quiver(0,0,0,0,0,1.28,color=navy,arrow_length_ratio=.08,linewidth=1.6)
ax3.text(0,0,1.36,'表面法向',ha='center',color=navy)
ax3.scatter([0],[0],[0],color=navy,s=18)
xx,yy=np.meshgrid(np.linspace(-.6,.6,2),np.linspace(-.6,.6,2))
ax3.plot_surface(xx,yy,np.zeros_like(xx),color=grey,alpha=.35,shade=False)
ax3.set(xlim=(-.65,.65),ylim=(-.65,.65),zlim=(0,1.45))
ax3.set_box_aspect((1,1,1.2))
ax3.view_init(elev=17,azim=-65)
ax3.set_axis_off()
fig.text(.55,.93,'(b) 角锥：法向周围的三维方向集合',fontsize=14,color=navy)
fig.text(.55,.08,'15°、30°均指半角；覆盖全部方位角 φ，不只是一个截面。',fontsize=10,color='#556877')
fig.text(.55,.035,'此图解释评价区域，不表示已经算得的发射图样。',fontsize=10,color='#556877')
fig.savefig(OUT/'geometry_and_cones.png',dpi=200)
fig.savefig(OUT/'geometry_and_cones.pdf')
print(OUT/'geometry_and_cones.png')
