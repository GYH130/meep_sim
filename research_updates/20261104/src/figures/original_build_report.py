from pathlib import Path
import json, hashlib
from datetime import datetime
from docx import Document
from docx.shared import Cm, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT = Path('<WORKSPACE_ROOT>')
OUT = ROOT/'outputs/teacher_report_20261104'
OUT.mkdir(parents=True, exist_ok=True)
FIG = OUT/'figures'
SUMMARY = ROOT/'outputs/ml_remaining4_mem24r2_20261103/all9_completed_summary_20261104.json'
S = json.loads(SUMMARY.read_text())
ROWS = S['individual_results']
M = [r['metrics'] for r in ROWS]
MLDIR = ROOT/'outputs/ml_surrogate_pilot_20261102_update_20261104'
ML = json.loads((MLDIR/'update_summary.json').read_text())
DA = json.loads((MLDIR/'data_audit.json').read_text())
PC = json.loads((MLDIR/'prospective_prediction_validation.json').read_text())
doc = Document()
sec = doc.sections[0]
sec.page_width,sec.page_height=Cm(21),Cm(29.7)
sec.top_margin,sec.bottom_margin=Cm(1.65),Cm(1.6)
sec.left_margin,sec.right_margin=Cm(2),Cm(2)
sec.header_distance,sec.footer_distance=Cm(.7),Cm(.7)

def font(style, size, bold=False):
    st = doc.styles[style]
    st.font.name='Arial Unicode MS'
    st.font.size=Pt(size)
    st.font.bold=bold
    st.font.italic=False
    st.font.color.rgb=RGBColor.from_string('000000')
    rf=st.element.get_or_add_rPr().rFonts
    for key in list(rf.attrib):
        if 'theme' in key.lower():del rf.attrib[key]
    rf.set(qn('w:eastAsia'),'Arial Unicode MS')
    st.paragraph_format.space_after=Pt(6)
    st.paragraph_format.line_spacing=1.22
for name,size,b in [('Normal',11,False),('Title',23,True),('Subtitle',11,False),('Heading 1',16,True),('Heading 2',12,True),('Caption',9,False)]:font(name,size,b)
for st in doc.styles:
    for border in st.element.xpath('.//w:pBdr'):
        border.getparent().remove(border)
doc.styles['Normal'].paragraph_format.widow_control=True
doc.styles['Heading 1'].paragraph_format.space_before=Pt(0)
doc.styles['Heading 2'].paragraph_format.space_before=Pt(7)
for name in ['Title','Heading 1','Heading 2']:
    doc.styles[name].paragraph_format.keep_with_next=True

h=sec.header.paragraphs[0]
h.text='钛基微腔定向热发射研究  ·  几何与响应图解'
h.runs[0].font.size=Pt(8)
h.runs[0].font.color.rgb=RGBColor.from_string('000000')
f=sec.footer.paragraphs[0];f.alignment=WD_ALIGN_PARAGRAPH.RIGHT
f.add_run('2026年11月4日  |  探索性结果  |  ')
fld=OxmlElement('w:fldSimple');fld.set(qn('w:instr'),'PAGE');f._p.append(fld)
for r in f.runs:r.font.size=Pt(8)

def p(t='',bold=False,style=None):
    z=doc.add_paragraph(t,style)
    if bold:
        for r in z.runs:r.bold=True
    return z
def heading(t):doc.add_heading(t,1)
def sub(t):doc.add_heading(t,2)
def page():doc.add_page_break()
def fig(name,caption,width=16.5):
    z=doc.add_paragraph();z.alignment=WD_ALIGN_PARAGRAPH.CENTER
    z.paragraph_format.space_after=Pt(2)
    z.add_run().add_picture(str(Path(name) if Path(name).is_absolute() else FIG/name),width=Cm(width))
    c=p(caption,style='Caption');c.paragraph_format.keep_with_next=False
def table(headers,rows,widths=None,size=9):
    t=doc.add_table(rows=1,cols=len(headers));t.autofit=False
    t.alignment=1
    if widths:
        for col,w in zip(t.columns,widths):col.width=Cm(w)
    pr=t._tbl.tblPr
    borders=OxmlElement('w:tblBorders')
    for k in ['top','left','bottom','right','insideH','insideV']:
        z=OxmlElement('w:'+k);z.set(qn('w:val'),'single');z.set(qn('w:sz'),'4');z.set(qn('w:color'),'D9D9D9');borders.append(z)
    pr.append(borders)
    mar=OxmlElement('w:tblCellMar')
    for k in ['top','bottom','left','right']:
        z=OxmlElement('w:'+k);z.set(qn('w:w'),'75' if k in ['top','bottom'] else '95');z.set(qn('w:type'),'dxa');mar.append(z)
    pr.append(mar)
    for i,v in enumerate(headers):t.rows[0].cells[i].text=str(v)
    repeat=OxmlElement('w:tblHeader');t.rows[0]._tr.get_or_add_trPr().append(repeat)
    for row in rows:
        cells=t.add_row().cells
        for i,v in enumerate(row):cells[i].text=str(v)
    for ri,row in enumerate(t.rows):
        trpr=row._tr.get_or_add_trPr();ns=OxmlElement('w:cantSplit');trpr.append(ns)
        for ci,c in enumerate(row.cells):
            if widths:c.width=Cm(widths[ci])
            tcpr=c._tc.get_or_add_tcPr();shade=OxmlElement('w:shd')
            shade.set(qn('w:fill'),'263F50' if ri==0 else ('F2F6F8' if ri%2 else 'FFFFFF'));tcpr.append(shade)
            for z in c.paragraphs:
                z.paragraph_format.space_after=Pt(0);z.paragraph_format.line_spacing=1.12
                for r in z.runs:
                    r.font.size=Pt(size);r.font.bold=ri==0
                    r.font.color.rgb=RGBColor.from_string('FFFFFF' if ri==0 else '111111')
    z=p();z.paragraph_format.space_after=Pt(0);z.paragraph_format.line_spacing=Pt(3);z.paragraph_format.space_before=Pt(0)
    return t

# Native editable Word equations. These use OMML fractions, integrals and scripts.
def mr(t):
    z=OxmlElement('m:r');pr=OxmlElement('w:rPr')
    fs=OxmlElement('w:sz');fs.set(qn('w:val'),'22');pr.append(fs);z.append(pr)
    a=OxmlElement('m:t');a.text=t;z.append(a);return z
def parts(parent, seq):
    for x in seq:
        parent.append(mr(x) if isinstance(x,str) else x)
    return parent
def subscript(base, sub):
    z=OxmlElement('m:sSub')
    z.append(parts(OxmlElement('m:e'),[base]));z.append(parts(OxmlElement('m:sub'),[sub]));return z
def superscript(base, sup):
    z=OxmlElement('m:sSup');z.append(parts(OxmlElement('m:e'),[base]));z.append(parts(OxmlElement('m:sup'),[sup]));return z
def frac(n,d):
    z=OxmlElement('m:f');z.append(parts(OxmlElement('m:num'),n));z.append(parts(OxmlElement('m:den'),d));return z
def mathbox(seq):
    z=OxmlElement('m:box');z.append(parts(OxmlElement('m:e'),seq));return z
def integ(lo,hi,body):
    z=OxmlElement('m:nary');pr=OxmlElement('m:naryPr')
    ch=OxmlElement('m:chr');ch.set(qn('m:val'),'∫');pr.append(ch)
    lc=OxmlElement('m:limLoc');lc.set(qn('m:val'),'subSup');pr.append(lc);z.append(pr)
    for name,value in [('supHide',hi),('subHide',lo)]:
        if not value:
            hide=OxmlElement('m:'+name);hide.set(qn('m:val'),'1');pr.append(hide)
    z.append(parts(OxmlElement('m:sub'),[lo]));z.append(parts(OxmlElement('m:sup'),[hi]));z.append(parts(OxmlElement('m:e'),body));return z
def eq(items,number):
    z=p();z.alignment=WD_ALIGN_PARAGRAPH.CENTER
    z.paragraph_format.space_before=Pt(3);z.paragraph_format.space_after=Pt(5)
    z.paragraph_format.line_spacing=1.0
    math=OxmlElement('m:oMath');parts(math,items);z._p.append(math)
    z.add_run(f'    ({number})').font.size=Pt(9)
    return z

V2=ROOT/'outputs/ml_batch_report_20261104_v2/figures'
FAB=ROOT/'outputs/fabrication_research_20261104/figures'
h.text='TA2钛基微结构定向热发射  ·  阶段研究报告'
def note(t):return p(t,style='Caption')

doc.add_paragraph('TA2钛基微结构\n定向热发射阶段研究报告','Title')
p('从二维机理预筛到三维微腔、膜层及机器学习辅助选点',style='Subtitle')
p('数据截止：2026年11月4日｜阶段总结与后续研究方案',style='Subtitle')
p('研究目标：面向8–13 μm大气窗口，探索激光可加工钛基结构的特征尺寸、介质填充与方向发射之间的关系，兼顾带内发射能力、目标角域功率和实际加工可实现性。',True)
fig('research_route.png','图1  研究路径与证据层级。二维全波、射线预筛、三维有限元和代理模型分别回答不同问题；浅色部分为后续任务。',17)
sub('本阶段最重要的进展')
p('已建立可追溯的“结构—材料—方向/偏振—功率审查—数据建模”工作链。三维底填充研究归档有91条求解记录、88个去重条件，覆盖8组独立高度–厚度组合；最新九项补样及保存结果功率审查已完成。')
p('主要发现：二维占空比存在局部优选；三维底填充厚度可显著增强单点吸收；增加腔高同时提高法向与斜向吸收，未显示更强法向集中。代理模型已可辅助补点，但设计覆盖仍需扩充。')
p('目前支持局部尺寸与填充效应比较、稀疏方向响应判断和代理选点；尚未完成8–13 μm带内及完整三维角锥性能验证。报告保留诊断与失败记录，将物理解算完成、方法审查和正式合格分开。')

page();heading('1  目标、评价量与方向性')
fig(FAB/'reciprocal_angle.png','图2  用吸收求解评估方向热发射。入射与出射通道互易；实际热发射实验不需要外部照明。图为原理示意，不是计算得到的发射瓣。',16.5)
p('对线性、无源、互易且满足适当热平衡条件的结构，可在相同波长及对应偏振通道下用吸收率求方向发射率。斜入射是计算方向响应的手段；最终目标仍是出射角域性能。')
eq([subscript('ε','out'),'(λ,−',subscript('k','in'),',',subscript('q','rev'),') = ',subscript('A','in'),'(λ,',subscript('k','in'),',q)'],1)
eq(['ε̄ = ½(',subscript('ε','s'),' + ',subscript('ε','p'),'),    ',subscript('k','in'),'/k = (sinθ cosφ, sinθ sinφ, −cosθ)'],2)
p('θ相对表面法线，φ为方位角。互易出射方向是−kin，方位角相差180°；理想A3具有六重对称，180°旋转等价，不对称斜槽不能直接沿用此等价。法向s/p表示代码选定的两正交线偏振。')
eq([subscript('ε','8–13'),'(θ,φ;T) = ',frac([integ('8 μm','13 μm',['ε̄(λ,θ,φ) ',subscript('B','λ'),'(T) dλ'])],[integ('8 μm','13 μm',[subscript('B','λ'),'(T) dλ'])])],3)
eq([subscript('B','λ'),'(T) = ',frac(['2 h c²'],[superscript('λ','5'),' (exp(',frac(['h c'],['λ ',subscript('k','B'),' T']),') − 1)'])],'3a')
p('Bλ(T)为普朗克谱辐亮度。波段指标必须指定温度；带内平均发射率与10.5 μm单点吸收率不同。当前重点是补齐可信结构响应，再计算波段及角域积分。')

page();heading('2  从离散方向到完整角域')
fig(V2/'cone_definition.png','图3  法向目标角锥与离散采样。右图示例仅指最新九项的两个方向，历史角点见图7。15°/30°为可调整评价半角，并非实测波束宽度。',14)
eq([subscript('P','Ω'),' = ',subscript('A','c'),' ',integ('8','13',[integ('Ω','',['ε̄ ',subscript('B','λ'),'(T) cosθ dΩ dλ'])]),',    dΩ = sinθ dθ dφ'],4)
eq([subscript('η','cone'),' = ',frac([subscript('P','cone')],[subscript('P','hemisphere')]),',    ',subscript('Q','cone'),' = ',frac([subscript('P','cone')],[subscript('P','blackbody,hemisphere')])],5)
p('Ac为投影单元面积。ηcone衡量本结构功率有多少落在目标锥内；Qcone衡量锥内功率相对同投影面积黑体半球功率的大小。二者同时报告，避免将总体发射很弱但角度比值较高的结构误判为优胜者。式(4)波长上下限为μm，数值计算须统一SI单位。')
table(['评价层级','本阶段状态','后续取得方式'],[
 ['单点/偏振吸收','已有真实三维结果','统一源、模式与功率定义'],
 ['角分辨响应','已有若干θ、φ离散点','同尺寸补s/p及关键方位角'],
 ['完整三维角锥功率','尚无全波积分结论','θ–φ覆盖及角度加密检验'],
 ['8–13 μm带内方向性','尚未完成','同结构的波长–角域联合积分']],[4.3,5.2,7.5],9.3)
p('先用稀疏角点识别法向集中或离轴峰值，再在峰值、半高处及角锥边缘加密。数据不足时展示真实散点，不补画连续发射瓣。')

page();heading('3  方法选择与已完成工作')
table(['方法','用途与已完成内容','适用边界'],[
 ['Meep / FDTD','二维斜槽尺寸及偏振预筛；平面基准、停止条件与功率诊断；三维平面可行性','保留波动效应；大尺度三维腔体细网格成本高'],
 ['几何光学射线','六类结构、尺寸与角域权衡预筛','忽略衍射/干涉/近场；灰体活性面是设计假设'],
 ['COMSOL / FEM','A3三维腔体、钛阻抗边界、SiO₂底填充及多方向/波长求解','非均匀网格适配几何；内存、模式与功率定义需核对'],
 ['TMM / 平面基准','平面膜层与材料输入交叉检查','不能替代三维腔体'],
 ['GP / 岭回归','质量标签、整结构留出验证、下一采样点建议','小样本代理，不能替代全波求解']],[3,8.2,5.8],9.3)
sub('为何由FDTD推进到FEM')
p('几十至上百微米腔体与中红外波长共存，使三维FDTD均匀细网格及时间步推进成本很高。FEM可沿几何和膜层局部布网，适合当前有限波长的频域微腔任务；并非此类结构在原理上只能用COMSOL。')
table(['数据集','规模','计数口径'],[
 ['三维底填充研究归档','91原始 / 88去重条件','8组H–t设计；不含所有控制/失败尝试'],
 ['其中保存场方法审查','37通过 / 14未通过 / 37待审查','按88个独立条件；正式合格数为0'],
 ['当前ML数据池','89原始 / 86去重条件','另2个9 μm警告点单列，未静默进入训练'],
 ['最新补样','9项求解及功率方法审查完成','增加交叉尺寸及偏振配对覆盖']],[5.1,5,6.9],9.1)
p('射线调用、重复网格、失败重试及模型预测分别计数，不充作新增物理结构。归档保留Meep、平面控制与所有重要诊断阶段。')

page();heading('4  Meep二维：占空比与离轴响应')
fig('meep_key_relationships.png','图4  真实二维斜槽及其响应：固定槽轴长和倾角，比较占空比；固定F=0.65比较三个波长。散点为已求解结果，连线不代表连续最优曲线。',17)
p('二维结构采用周期50 μm、槽轴长40 μm、倾角30°、钛厚60 μm。F=0.50、0.65、0.80对应开口25、32.5、40 μm；槽深与剩余钛厚也随之改变，不能将占空比效应理解为仅改变一项局部宽度。')
p('10.5 μm、p偏振下，+30°净通量吸收率依次约13.33%、19.97%、16.10%；其+30°与−30°差约4.82、11.04、7.82个百分点。图中主曲线采用独立体吸收，两种方法给出相同局部优选：F=0.65是三个样本中较好候选，不是全局最优。')
p('F=0.65时，8、10.5、13 μm的±30°差分别约11.42、11.04、8.56个百分点，支持继续研究非对称槽的方向差。结论限于二维、一个偏振及少量波长，不能直接推广到三维非偏振带内角域。')
note('当前结构10点含8条探索收敛、2条基础结果；12个配套平面控制中4通过、8存在反射精度问题。结构通量/体吸收偏差约0.269–0.380个百分点；方向差的方法差约0.00773个百分点，不是总误差。')

page();heading('5  六结构及尺寸预筛')
fig('structure_screen.png','图5  六类几何及射线预筛权衡。灰体活性面假设与真实材料全波模型分开；Q₁₅和C₁₅分别表征锥内绝对功率及锥内占比。',17)
p('初筛比较曲壁腔A0、底部偏移腔A1/A2、直壁腔A3及B/C不对称结构，再筛深度和底部尺度。A系列射线基准采用周期80 μm、高90 μm、开口对边69.5456 μm、底部20 μm；与后续FEM底部15 μm并非同一几何。')
p('M2模型为钛壁面配合假定吸收率0.9的活性面，它不是已实测SiO₂。密集候选中15°锥内占比约18–22%，包含50%半球功率的锥半角约27–35°；该积分半角不是角分布半高宽。')
p('增大深度可以增加目标角内绝对功率，也可能增强锥外壁面辐射、降低锥内占比。曲壁与直壁候选存在不同权衡；选择A3进入全波验证同时考虑激光加工及建模可控性。')
note('六结构主筛批次有114,180次调用、31,140个波段–角度样本，不是同等数量的独立实验或全波算例。射线预筛不验证衍射、干涉或波长级定向共振。')

page();heading('6  三维全波模型与材料输入')
fig(V2/'geometry_atlas.png','图6  A3六角截锥腔：三维剖切、周期俯视与垂直截面。示例H100/t2，SiO₂仅底部填充，侧壁不覆膜。',17)
table(['量或材料','当前模型','物理含义'],[
 ['周期 / 开口 / 底部对边距','80 / 69.5456 / 15 μm','横向尺寸固定，六重对称'],
 ['腔高H / 填充t','已采样H=60–120、t=1–3 μm子集','改变H同时改变壁面斜率及面积'],
 ['钛','可追溯光学钛常数 + 阻抗边界','不是用户TA2试样实测常数'],
 ['SiO₂','Kischkat等薄膜n,k；显式填充域','原数据为约510 nm溅射膜，不等于多孔溶胶或颗粒膜']],[4.8,6,6.2],9.2)
p('模型采用周期边界与入射/反射模式功率核算。钛阻抗边界近似不透射；SiO₂保留体损耗。平滑壁面、尖角及均匀填充均为理想化假设，后续须用实际形貌和膜层光学常数替换。')

page();heading('7  历史三维结果：已有角点，尚未闭合角域')
fig('historical_comsol_scope.png','图7  受控角度/波长点与ML设计覆盖。●为功率方法已审核；×为该比较未通过复核；空心谱点待审查。覆盖面板是86条件ML池，另2项9 μm警告结果单列。',16)
p('H60/t3与H100/t3在10.5 μm、30°及60°已有s/p对。H100非偏振吸收约18.34%、12.25%，高于H60约14.04%、8.63%；13 μm法向约41.81%、38.80%。这些成对结果通过保存场功率方法检查。')
p('部分法向10.5 μm及斜向13 μm点有方法复核问题，保留作敏感性线索，不与通过点等权拟合规律。H80/t3法向p在9、10.5、11.5、13 μm有四个谱点；仍不是收敛连续光谱或带内积分。')
p('历史已包含θ=0°、15°、22.5°、30°、60°及部分φ样本，说明三维角分辨模拟已开展。缺口是同结构、同质量标准的完整s/p与θ–φ–λ覆盖，而不是软件不能模拟方向性。')
p('旧H120法向10.5 μm非偏振总吸收在t1/t2/t3约36.55%/31.37%/37.06%，其方法质量尚未统一；13 μm法向p从t1到t3约24.77%→42.74%。平面TMM的8–13 μm未加权均值在t1.95/t2约54.814%/54.809%（不是腔体）。这些线索说明膜厚效应依赖波长及几何，不能沿用一个通用最优厚度。')

page();heading('8  最新九项：局部尺寸与填充效应')
fig(V2/'response_relationships.png','图8  最新九项10.5 μm结果：高度响应、填充吸收分解、偏振配对及两点方向比。细线只连接实际端点，没有外推。',17)
p('固定t2，H80→H100使法向非偏振吸收由26.86%增至29.52%，30°/30°由15.22%增至17.92%，提高2.655与2.705个百分点。A30/A0由0.5665升至0.6072，斜向相对法向也变强，不能据此宣称法向更集中。')
p('固定H80、法向s，t2→t3使吸收由26.87%升至34.13%，提高7.264个百分点，约94%的增量来自SiO₂积分吸收。两个端点支持继续探索填充厚度，但不能证明越厚越好。')
p('法向s/p差约0.003–0.005个百分点；斜向差达1.603与3.283个百分点。非偏振角域评价应保留s/p配对，不能只选择高响应偏振。')

page();heading('9  功率核算、误差与可信度')
eq(['R = ',frac(['Σ ',superscript(subscript('P','m,q'),'out')],[subscript('P','in')]),',    ',subscript('A','port'),' = 1 − R − T,    T ≈ 0'],6)
eq([subscript('A','loss'),' = ',frac([integ('Vfilm','',['Q dV']),' + ',integ('STi','',['½ Re(',subscript('Z','s'),') |',subscript('H','t'),'|² dS'])],[subscript('P','in')])],7)
eq(['Q = ½ ω ',subscript('ε','0'),' ',subscript('ε','loss'),'″ |E|²,    ',subscript('ε','loss'),'″ = 2nk ≥ 0'],8)
fig(V2/'numerical_agreement.png','图9  最新九项多种功率定义及保存结果独立核算。材料损耗、端口差与原始净通量分别保留；残差不是完整总不确定度。',16)
p('九项保存结果功率方法均通过，最大绝对残差0.173个百分点、相对0.829%。原始A_flux仍比A_loss高约1.07–2.76个百分点，原始警告没有改为“通过”。求解完成、保存场核算通过和正式网格合格须分开。')
p('原0.002及探索0.005绝对容差、2%相对核算检查均保留。探索条件不是总误差足够小的证明。提高吸收会减小同等绝对误差的相对值，但不能消除源、模式、功率定义或几何错误。')
note('式(8)使用正损耗量εloss″，复介电常数符号须与求解器时间约定一致；R覆盖传播级次及两偏振。最新九项仍无正式网格收敛证明。')

page();heading('10  加工与填充：从理想模型到试样')
fig(FAB/'fill_morphologies.png','图10  裸腔、理想底填充、共形覆膜及实际干燥沉积的区别。图为工艺形貌示意；不同形貌不能共用未经测量的厚度与n,k。',17)
p('文献证明TA2/钛合金的上百微米阵列和深腔具有激光加工可行性，但没有证明本设备能直接得到当前理想尖角、光滑壁面及窄底密排腔。H80/H100侧壁相对竖直约18.8°/15.3°，顶面间隙约10.45 μm；再沉积、氧化、圆角和粗糙度可能改变响应。')
p('SiO₂是可追溯、具有中红外声子损耗的对照材料，尚非已证明最优。溶胶–凝胶/局部滴涂便于低成本试样，但会产生收缩、裂纹、咖啡环及孔隙；ALD趋向共形并覆盖侧壁；PVD有视线遮挡。几种路线中的“2 μm”并不等价。')
p('先做裸TA2宽松阵列/孤立腔和平面见证片，取得实际截面、氧化及粗糙度，再对比薄覆膜和局部填充。用FTIR及合适的材料表征更新模型。H100/t2理想干态体积仅约0.419 pL/腔，滴量不能直接当作干膜体积。')
note('工艺依据涵盖TA2蜂窝激光加工、钛箔漏斗孔和钛合金微铣削；本设备配方尚未验证。Kischkat溅射膜n,k不能直接验证多孔溶胶SiO₂。')

page();heading('11  机器学习：已完成验证及当前定位')
fig(V2/'ml_evidence.png','图11  新H100/t2结构的事先预测与实际计算，以及H–t覆盖。四个方向/偏振条件仍只对应一组独立几何。',17)
p('训练已在Mac的Python环境实际执行，使用NumPy、SciPy、scikit-learn及joblib，无需独立图形软件。Windows服务器上的COMSOL提供标签。模型输入x=(H,t)，按波长、方向和偏振分条件拟合材料总吸收。')
eq([subscript('μ','*'),' = ',subscript('μ','y'),' + ',superscript(subscript('k','*'),'T'),' ',superscript(mathbox(['(K + ',superscript(subscript('σ','n'),'2'),' I)']),'−1'),' (y − ',subscript('μ','y'),')'],9)
eq([superscript('MAE','pp'),' = ',frac(['100'],['N']),' Σ |',subscript('Â','i'),' − ',subscript('A','i'),'|'],10)
p('GP使用Matérn 5/2核，保留均值基线及岭回归。H100/t2事先预测误差依次为2.54、1.64、1.12个百分点，提示局部选点价值；验证对象只有一个新结构。历史整结构留出中法向GP约4.20个百分点，未稳定优于均值约4.12。')
p('当前8组设计仍稀疏；26个10.5 μm条件含13条已审核、9条待审核、4条方法失败诊断。下一步先交叉补样、保留质量标签，以整套H/t留出防止数据泄漏。当前定位是代理预测与主动选点，不是直接端到端逆向设计。')

page();heading('12  后续方案、优先级与工作量')
table(['优先级 / 任务','目的及交付','判断依据 / 工作量'],[
 ['P1 可加工裸腔及填充试样','截面、壁面/氧化、平面见证片与FTIR；更新几何及n,k','优先消除理想底填充与实际形貌差；周期依设备/工艺'],
 ['P1 H60/t2、H80/t1补样','相同材料/横向几何；10.5 μm法向、30°/30°的s/p','增加2组设计、8条件；暂估24–48小时'],
 ['P2 优选两结构角域','先稀疏θ/φ，补峰值、边缘及s/p；计算目标锥指标','先单波长方向图，再检验积分稳定性'],
 ['P2 8–13 μm及后续3–20 μm','同结构补光谱，计算带内及锥内功率','只对可信候选扩谱；材料范围不等于已求解范围'],
 ['P3 代理主动学习','整结构留出，对比基线，按预测与不确定性选点','减少求解次数；最终候选仍需全波及实验验证']],[4.2,7.5,5.3],9.3)
sub('计算组织与预计时间')
p('本批九项已结束，下一批尚未启动。按本批单项约3–6小时估计，两项小批约6–12小时，两个新结构完整8项约24–48小时。实际网格、材料版本及共享内存会改变耗时，目前没有新的北京时间完工窗口。')
p('默认串行32核；启动可用RAM≥100 GiB、运行保留64 GiB，本项目RSS软112/硬128 GiB；不改动他人ORCA/FDTD任务。优先独立诊断和局部修订，避免为相同残差重跑全部求解，保留失败原因和版本指纹。')
sub('进入大规模寻优前的里程碑')
p('建立一条实际可加工路线、一个独立功率核对且具有局部网格敏感性记录的材料/形貌模型，以及同结构波长–方向–偏振采样。之后让代理模型在加工约束内围绕带内发射能力、目标角域功率及工艺难度多目标选点。')

page();heading('附表A  最新九项核心数值')
table(['编号','H/t\nμm','θ/φ','偏振','R\n%','A_loss\n%','SiO₂\n%','Ti\n%'],[
 [f'ML{i}',f"{m['height_um']:.0f}/{m['floor_fill_thickness_um']:.0f}",'0°' if m['theta_deg']==0 else '30°/30°',m['polarization'],f"{m['R_total']*100:.3f}",f"{m['A_loss_total']*100:.3f}",f"{m['A_film']*100:.3f}",f"{m['A_Ti_H']*100:.3f}"] for i,m in enumerate(M,1)
],[1.3,1.7,2.3,1.2,2.5,2.7,2.6,2.7],8.8)
p('A_loss=SiO₂+Ti；A_port=1−R。各项λ=10.5 μm，均已完成保存结果功率方法审查；原始净通量警告仍保留，没有正式网格合格声明。队列于北京时间11月4日13:16:35结束。')
sub('新增遗漏记录已补入研究台账')
p('H100/t3、9 μm法向p/s的A_loss约38.852%/38.870%，两项于11月1日完成，但方法审核失败，另表保留。研究归档91/88与ML池89/86因此不同；未把这些诊断点静默用于模型训练。')
sub('其他阶段的结果边界')
p('旧二维覆膜系列0/1/3 μm曾得约7.03%/50.24%/81.84%均值，但仅是9、10、12 μm三点在9–12 μm的离散积分，不是8–13 μm连续性能，也不是A3三维腔体。Meep三维平面r32测试吸收约3.727%，对应解析约3.798%；没有同腔体的三维Meep验证。')
sub('资源与后台任务记录')
p('本批首三项峰值RSS47.98 GiB，最后四项24 GiB OOC缓存队列峰值34.09 GiB。11月4日14:26:57快照：本项目无活动进程、可用内存180.67 GiB、CPU约18%、D盘健康。此为带时间戳记录，不是阅读时的实时负载。')

page();heading('附录B  数据、代码与主要依据')
p('公开归档：GYH130/meep_sim，codex/research-summary-20261104分支的research_updates/20261104目录。包含本报告、图表源程序/矢量图、汇总数据、求解与审查核心实现、代理模型代码。逐文件SHA256及实际Git提交记录随发布清单保存，历史分支保留。')
p('公开包不含SSH/网络凭据、主机用户信息、虚拟环境、大型MPH或全量三维场。原始场本地保留。数据表保留原始accepted、保存场审查、警告及正式资格，读者可区分求解、审核和结论。')
p('图形借鉴经典论文的“真实几何/剖面—尺寸材料—参数小多图—一致性/权衡”组织；统一面板字母、颜色、单位与偏振，保留原始散点，输出SVG及高分辨率PNG。公式为可编辑Word公式。不复制他人成果，不用插值补齐未算区域。')
refs=[
 ('图形组织','Raman et al., Nature 515, 540–544 (2014), Fig.1–3.','10.1038/nature13883'),
 ('结构–角度–偏振图组','Greffet et al., Nature 416, 61–64 (2002), Fig.1–4.','10.1038/416061a'),
 ('方向发射及互易','Miller, Zhu & Fan, PNAS 114, 4336–4341 (2017).','10.1073/pnas.1701606114'),
 ('SiO₂光学数据','Kischkat et al., Applied Optics 51, 6789–6798 (2012).','10.1364/AO.51.006789'),
 ('学习图形组织','Peurifoy et al., Science Advances 4, eaar4206 (2018), Fig.1–2.','10.1126/sciadv.aar4206'),
 ('TA2激光蜂窝','Wang et al., Optik (2019), article 162953.','10.1016/j.ijleo.2019.162953'),
 ('钛激光漏斗孔','Ye et al., Scientific Reports 6, 37591 (2016).','10.1038/srep37591'),
 ('钛合金微铣削','Calabrese et al., JMMP 5, 125 (2021).','10.3390/jmmp5040125'),
 ('干燥沉积形貌','Shi et al., ACS Applied Materials & Interfaces (2019).','10.1021/acsami.8b21731'),
 ('SiO₂共形ALD','Dingemans et al., Journal of the Electrochemical Society (2012).','10.1149/2.067203jes')]
for topic,title,doi in refs:
 z=p(f'{topic}：{title}  https://doi.org/{doi}');z.paragraph_format.line_spacing=1.05;z.paragraph_format.space_after=Pt(5)
 for r in z.runs:r.font.size=Pt(8.5)
note('文献工艺可行性不等于本设备、本形貌已实现。详细材料/工艺证据及源索引保存在归档数据中。')

FINAL=OUT/'final';FINAL.mkdir(exist_ok=True)
TARGET=FINAL/'TA2钛基微结构定向热发射_阶段研究报告_20261104.docx'
doc.save(TARGET)
manifest={'report':str(TARGET.relative_to(ROOT)),'cutoff':'2026-11-04','physical_simulations_started_this_report':0,'all9_summary_sha256':hashlib.sha256(SUMMARY.read_bytes()).hexdigest(),'native_equation_count':len(doc.element.xpath('.//m:oMath')),'figure_count':len(doc.inline_shapes),'status':'built_pending_render_qa'}
(OUT/'report_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
print(json.dumps(manifest,ensure_ascii=False))

