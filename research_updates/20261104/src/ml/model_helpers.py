"""Exploratory condition-specific surrogate pilot. Never launches a solver.

Every held-out design is removed from all condition-specific training sets.
Quality flags and alternative power labels remain separate observations.
"""
from pathlib import Path
from collections import defaultdict
import datetime as dt
import hashlib
import html
import json
import time
import warnings
import platform
import numpy as np
import scipy
from scipy.spatial import Delaunay
from scipy.stats import spearmanr
import sklearn
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
SEED = 20261102
LABELS = ['A_loss_total', 'A_port']
CANDIDATES = [(80., 2.), (100., 2.), (60., 2.), (80., 1.)]

def write(name, obj):
    (HERE / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def cond(row):
    return f"l{row['wavelength_um']:g}_th{row['theta_deg']:g}_phi{row['phi_deg']:g}_{row['polarization']}"

def features(rows):
    return np.array([[r['height_um'], r['film_um']] for r in rows], dtype=float)

def fit(rows, label, method):
    x = features(rows)
    y = np.array([r[label] for r in rows], dtype=float)
    scaler = StandardScaler().fit(x)
    z = scaler.transform(x)
    if method == 'mean':
        estimator = DummyRegressor(strategy='mean')
    elif method == 'linear_ridge':
        # Fixed modest regularization: not selected using held-out predictions.
        estimator = Ridge(alpha=0.1)
    else:
        # Isotropic kernel limits hyperparameter count in this tiny dataset.
        kernel = ConstantKernel(1., (1e-3, 1e3)) * Matern(1., (0.05, 20.), nu=2.5)
        kernel += WhiteKernel(0.01, (1e-5, 1.))
        estimator = GaussianProcessRegressor(kernel=kernel, normalize_y=True,
            alpha=1e-10, random_state=SEED, n_restarts_optimizer=2)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        estimator.fit(z, y)
    return {'scaler': scaler, 'estimator': estimator,
            'warnings': sorted(set(str(w.message) for w in caught)),
            'train_designs': sorted({r['design_id'] for r in rows}),
            'label': label, 'method': method,
            'status': 'EXPLORATORY_UNCALIBRATED_NOT_FORMAL'}

def predict(model, x):
    z = model['scaler'].transform(np.atleast_2d(x))
    if model['method'] == 'gp_matern':
        m, s = model['estimator'].predict(z, return_std=True)
        return m, s
    return model['estimator'].predict(z), None

def support(x, query):
    # Distance scaling is fixed from the proposed physical parameter range.
    z = (x - [60., 1.]) / [60., 2.]
    q = (np.array(query) - [60., 1.]) / [60., 2.]
    centered = z-z[0]
    rank = int(np.linalg.matrix_rank(centered))
    if rank == 2:
        inside = bool(Delaunay(z).find_simplex(q, tol=1e-10) >= 0)
    elif rank == 1:
        _, _, v = np.linalg.svd(centered, full_matrices=False)
        u = v[0]
        positions = centered@u
        distance = q-z[0]
        projection = float(distance@u)
        inside = bool(np.linalg.norm(distance-projection*u) < 1e-8 and
                      positions.min()-1e-8 <= projection <= positions.max()+1e-8)
    else:
        inside = bool(np.linalg.norm(q-z[0]) < 1e-8)
    return {'inside_training_convex_hull': inside, 'training_affine_rank': rank,
            'warning': None if inside else 'NEW_PARAMETER_REGION_REQUIRES_PHYSICAL_SOLVE'}

def main():
    started = time.monotonic()
    source = HERE / 'dataset.json'
    dataset = json.loads(source.read_text())
    raw = dataset['baseline_rows']
    rows = [r for r in raw if r['exploratory_fit_eligible']]
    assert len({(r['design_id'], cond(r)) for r in rows}) == len(rows)
    assert all(r['physical_solver_completed'] and not r['hard_failed_checks'] for r in rows)
    conditions = sorted({cond(r) for r in rows})
    groups = sorted({r['design_id'] for r in rows})
    methods = ['mean', 'linear_ridge', 'gp_matern']
    oof, metrics, candidates, checks, model_warnings = [], [], [], [], []
    (HERE/'models').mkdir(exist_ok=True)
    for label in LABELS:
        for condition in conditions:
            subset = [r for r in rows if cond(r) == condition and r[label] is not None
                      and np.isfinite(r[label])]
            for method in methods:
                for heldout in groups:
                    train = [r for r in subset if r['design_id'] != heldout]
                    test = [r for r in subset if r['design_id'] == heldout]
                    if not test:
                        continue
                    assert train and not ({r['design_id'] for r in train} & {heldout})
                    model = fit(train, label, method)
                    mean, std = predict(model, features(test))
                    model_warnings += [{'condition':condition,'label':label,'method':method,
                                        'heldout':heldout,'warnings':model['warnings']}] if model['warnings'] else []
                    for i, r in enumerate(test):
                        assert np.isfinite(mean[i])
                        item = {'condition':condition,'label':label,'method':method,
                                'heldout_design':heldout,'train_designs':model['train_designs'],
                                'case_id':r['case_id'],'observed':r[label],
                                'predicted':float(mean[i]),
                                'std':float(std[i]) if std is not None else None,
                                'quality_class':r['quality_class'],
                                'method_checked':r['method_checked_fit_eligible']}
                        oof.append(item)
                these = [p for p in oof if p['condition']==condition and p['label']==label and p['method']==method]
                y = np.array([p['observed'] for p in these]); p = np.array([p['predicted'] for p in these])
                e = p-y
                rho = float(spearmanr(y,p).statistic) if len(y)>=3 and np.ptp(y)>0 and np.ptp(p)>0 else None
                coverage = float(np.mean([abs(v['predicted']-v['observed']) <= 1.96*v['std'] for v in these])) if method=='gp_matern' else None
                metrics.append({'condition':condition,'label':label,'method':method,
                    'n_test_designs':len(these),'MAE_pp':float(np.mean(abs(e))*100),
                    'RMSE_pp':float(np.sqrt(np.mean(e**2))*100),
                    'max_abs_error_pp':float(max(abs(e))*100),'spearman':rho,
                    'nominal_95_interval_coverage':coverage,
                    'interval_warning':'Conditional GP interval; not material/systematic-error coverage or calibrated guarantee.'})
                final = fit(subset, label, method)
                final['condition'] = condition
                final['dataset_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
                joblib.dump(final, HERE/'models'/f'{condition}_{label}_{method}.joblib')
                if method == 'gp_matern':
                    for h,t in CANDIDATES:
                        m,s=predict(final, [[h,t]])
                        candidates.append({'design_id':f'H{h:g}_t{t:g}','height_um':h,'film_um':t,
                            'condition':condition,'label':label,'predicted':float(m[0]),
                            'conditional_std':float(s[0]),'interval95':[float(m[0]-1.96*s[0]),float(m[0]+1.96*s[0])],
                            'mean_outside_physical_bounds':bool(m[0]<0 or m[0]>1),
                            'selection_basis':'Coverage-driven proposed design, not validated optimum or BO acquisition.',
                            **support(features(subset),[h,t])})
            strict = [r for r in subset if r['method_checked_fit_eligible']]
            checks.append({'condition':condition,'label':label,'exploratory_count':len(subset),
                'method_checked_count':len(strict),
                'method_checked_only_fit':'NOT_FITTED: only 0-2 designs; no independent H/t validation possible',
                'excluded_from_exploratory_count':len([r for r in raw if cond(r)==condition])-len(subset)})
    write('heldout_predictions.json',oof)
    write('baseline_metrics.json',metrics)
    write('candidate_predictions.json',candidates)
    write('quality_sensitivity.json',{'subsets':checks,'model_warnings':model_warnings,
          'label_policy':'A_loss_total primary exploratory label; A_port separately fitted sensitivity channel; no mixed flux correction.',
          'no_formal_qualification_claim':True})
    # Coverage and genuine held-out errors, not a training-fit illustration.
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    designs=sorted({(r['height_um'],r['film_um']) for r in dataset['rows']})
    ax.scatter(*np.array(designs).T,s=100,label='Archived designs',color='#275c7b')
    ax.scatter([80,100],[2,2],s=130,facecolors='none',edgecolors='#c36a25',linewidths=2,label='Proposed first two designs')
    for h,t in designs: ax.annotate(f'H{h:g}/t{t:g}',(h,t),xytext=(5,8),textcoords='offset points',fontsize=9)
    ax.set(xlabel='Cavity height H (µm)',ylabel='Floorfill thickness t (µm)',xticks=[60,80,100,120],yticks=[1,2,3],xlim=(55,130),ylim=(.65,3.6),title='Six independent designs; sparse H × t coverage')
    ax.grid(alpha=.2); ax.legend(loc='lower left',fontsize=9)
    fig.tight_layout();fig.savefig(HERE/'design_coverage.png',dpi=170);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4.2))
    for ax,label in zip(axs,LABELS):
        for method,color in zip(methods,['#888888','#d28730','#275c7b']):
            z=[p for p in oof if p['label']==label and p['method']==method]
            ax.scatter([p['observed']*100 for p in z],[p['predicted']*100 for p in z],s=30,label=method,alpha=.8,color=color)
        lo,hi=0,60;ax.plot([lo,hi],[lo,hi],'k--',lw=1)
        ax.set(xlabel='Archived value (%)',ylabel='Held-out prediction (%)',title=label+' — leave one design out',xlim=(lo,hi),ylim=(lo,hi))
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Exploratory records with retained power-method warnings',fontsize=11)
    fig.tight_layout();fig.savefig(HERE/'heldout_validation.png',dpi=170);plt.close(fig)
    first=[]
    for h,t in CANDIDATES[:2]:
        for th,ph in [(0,0),(30,30)]:
            for pol in ['s','p']:
                first.append({'id':f'ML_H{h:g}_t{t:g}_l10p5_th{th}_phi{ph}_{pol}',
                    'height_um':h,'floorfill_um':t,'wavelength_um':10.5,'theta_deg':th,'phi_deg':ph,'polarization':pol})
    missing=[]
    for h in [80.,100.]:
        for th,ph in [(0.,0.),(30.,30.)]:
            for pol in ['s','p']:
                if not any(r['height_um']==h and r['film_um']==3 and r['theta_deg']==th and r['phi_deg']==ph and r['polarization']==pol for r in rows):
                    missing.append({'height_um':h,'floorfill_um':3.,'wavelength_um':10.5,'theta_deg':th,'phi_deg':ph,'polarization':pol})
    historical=[]
    root=HERE.parents[1]
    for rel in ['outputs/height_direction_20261028/completed_summary_20261029.json',
                'outputs/height_midangle_20261029/completed_summary_20261029.json']:
        for item in json.loads((root/rel).read_text())['individual_results']:
            historical.append({'case_id':item['case_id'],'wall_hours':item['case_wall_elapsed_s']/3600,'source':rel})
    first_block=missing+[c for c in first if c['height_um']==80 and c['theta_deg']==0]
    write('next_batch_plan.json',{'state':'PLANNED_NOT_DEPLOYED_NOT_LAUNCHED','primary_designs':first,
        'missing_matched_baselines':missing,'expected_new_physical_solves':len(first)+len(missing),
        'suggested_first_block':first_block,
        'first_block_purpose':'Complete H80/t3 normal s, then H80/t2 normal s/p; matched thickness comparison before angle extension.',
        'design_decision':'Controlled H/t cross comparisons; model predictions advisory only.',
        'first_step':'Reconcile existing queue, then prepare first H80/t2 case using established source/material and bounded resource controls.',
        'new_solver_time_limit':None,'new_cumulative_budget':None,
        'resources':'Reuse current authorized guards; measure live headroom and first-case footprint before allocation.',
        'historical_wall_time_evidence':historical,
        'scheduling_estimate':'Planning only: prior 10.5um cases span 0.87-3.29h wall time. Allow 5-8h for first 3, 12-30h for all 9; new mesh, export/review and shared resource wait can extend these ranges. No start time booked; re-estimate after first new case.',
        'scope':'A3 optical-Ti IBC plus SiO2 floorfill; not conformal film, measured TA2, full-band integral or cone efficiency.',
        'phase2':'For shortlisted designs, add spectral/angular samples for 8–13µm power and directional integrals.'})
    write('pilot_manifest.json',{'state':'LOCAL_EXPLORATORY_BASELINES_TRAINED','generated_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
        'elapsed_seconds':time.monotonic()-started,'python':platform.python_version(),
        'versions':{'numpy':np.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__,'matplotlib':matplotlib.__version__},
        'dataset_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'training_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'methods':methods,'conditions':conditions,'labels':LABELS,'n_baseline_rows':len(rows),
        'n_independent_designs':len(groups),'n_fitted_final_models':len(LABELS)*len(conditions)*len(methods),
        'new_physical_solves_started':0,'formal_qualification':False,
        'validation':'Leave one whole H/t design out across all conditions; all scalers fit on train only. No row-random CV. Fixed model choices; no test-set tuning.',
        'qa':{'no_duplicate_design_condition':True,'train_test_design_disjoint':True,'all_predictions_finite':True},
        'limits':['six designs only','uncertainty is model conditional and not calibrated','normal-condition quality warnings retained','no full 8–13µm cone labels','new region predictions require physical verification']})
    table=''.join('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in [m['condition'],m['method'],m['n_test_designs'],f"{m['MAE_pp']:.3f}",f"{m['max_abs_error_pp']:.3f}"] )+'</tr>' for m in metrics if m['label']=='A_loss_total')
    report=f'''<!doctype html><meta charset="utf-8"><title>代理模型首轮启动结果</title>
<style>body{{max-width:1050px;margin:40px auto;font:17px/1.7 system-ui;color:#20303d;padding:0 25px}}table{{border-collapse:collapse;width:100%;font-size:15px}}th,td{{padding:7px;border-bottom:1px solid #ddd;text-align:left}}img{{max-width:100%}}.note{{background:#fff4df;padding:18px}}a{{color:#276282}}</style>
<h1>小样本代理模型：首轮启动结果</h1>
<p>已完成本地数据整理与真实归档仿真数据训练。范围为 A3 光学钛阻抗边界、SiO₂底填充、10.5 μm。未部署或启动新物理解算。</p>
<p class="note">本轮是探索性基线。既有法向结果的功率核算警告、各算例正式资格与网格限制均保留。预测区间不是材料误差或数值系统误差的保证。</p>
<p>基线包含 {len(rows)} 个不同物理条件、{len(groups)} 组 H/t 设计。按完整结构留出验证；同一结构的其他条件不进入该折训练。比较均值、线性岭回归、高斯过程；另用 A_port 独立重做敏感性模型。</p>
<img src="design_coverage.png" alt="结构覆盖图"><h2>真正留出预测的误差</h2>
<p>表内误差为吸收率百分点。每个条件只有3至6个设计，不能据此宣称可靠宽域预测；按误差和实际补算结果判断是否值得采用GP。</p>
<p><b>本轮 GP 尚未展示预测收益。</b>法向 p/s 的 GP 留出平均误差分别为3.85/4.83个百分点，均值基线为3.56/4.36个百分点。30°样本只有3组且填充厚度均为3 μm，无法学习厚度变化。下一批先补交叉覆盖，模型作为辅助，不以它决定“最优结构”。</p>
<table><tr><th>条件</th><th>模型</th><th>留出设计数</th><th>MAE/百分点</th><th>最大误差/百分点</th></tr>{table}</table>
<img src="heldout_validation.png" alt="留出验证图"><h2>下一批建议</h2>
<p>先补 H80/t2 和 H100/t2，各做法向及 θ30°/φ30°的 s/p，共8项；补齐匹配 H80/t3 法向s基准则共9项。推荐依据是填补H/t交叉覆盖，不是模型证明它们最优。H80/t2位于现有法向样本凸包之外；所有t2斜向点也缺乏已有厚度变化支撑。</p>
<p>先安排3项：H80/t3法向s补基准，再做H80/t2法向s/p，形成同高度厚度对照；其余6项随后推进。已有10.5 μm算例实际用时约0.87–3.29小时/项，第一组暂估5–8小时，全部9项暂估12–30小时；新网格、导出审核和共享资源等待可能延长。未启动，因此没有北京时间完成时刻承诺。</p>
<p>A_port与A_loss_total来自同一物理解，是功率核算方法的敏感性对照，不是独立实验真值。GP白噪声项由数据拟合，不能当作COMSOL数值误差。补算与交叉验证后再决定是否使用主动学习选点。</p>
<p>并行整理实验“加工参数—实测尺寸/覆盖—匹配光谱/角谱”。第一版模型暂不预测激光加工参数，也不替代8–13 μm角域积分。</p>
<p><a href="baseline_metrics.json">全部验证指标</a> · <a href="quality_sensitivity.json">质量与模型警告</a> · <a href="candidate_predictions.json">候选预测及支持范围</a> · <a href="next_batch_plan.json">未启动的计算计划</a> · <a href="pilot_manifest.json">训练记录</a></p>'''
    (HERE/'startup_report.html').write_text(report)
    print(json.dumps({'state':'LOCAL_EXPLORATORY_BASELINES_TRAINED','rows':len(rows),'designs':len(groups),'models':24,'seconds':time.monotonic()-started,'new_solves':0}))

if __name__ == '__main__':
    main()
