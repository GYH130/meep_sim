"""Reproduce existing pilot and fit a pooled ray+FEM residual prototype.

Only local scalar data are read. FEM solvers and remote hosts are not called.
Whole H/b/t designs are held out; wavelengths/angles/polarizations are not designs.
"""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import csv, json, hashlib, time, datetime, platform, collections, resource
import numpy as np, scipy, sklearn, joblib
from surrogate import fit,predict,metrics,condition,METHODS

HERE=Path(__file__).resolve().parent
RAY=HERE/'inputs'
SOURCE=RAY/'dataset.json'

def dump(name,obj):
    (HERE/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def main():
    start=time.monotonic()
    (HERE/'inputs').mkdir(exist_ok=True);(HERE/'models').mkdir(exist_ok=True)
    evidence=[SOURCE,RAY/'ray_vs_fem.csv']
    evidence.extend(RAY/name for name in ['first_repaired_scalar_summary.json','second_repaired_scalar_summary.json'])
    evidence.extend(HERE/'ray_engine'/name for name in ['ray_material.cpp','raylib.py','materials.py','inputs/Ti_Ordal_material.json','inputs/SiO2_Kischkat.yml'])
    dataset=json.loads(SOURCE.read_text());old={r['case_id']:r for r in dataset['rows']}
    archived=list(csv.DictReader((RAY/'ray_vs_fem.csv').open()))
    rows=[]
    for a in archived:
        r=dict(old.get(a['case_id'],{}))
        for key in ['height_um','film_um','wavelength_um','theta_deg','phi_deg']:
            r[key]=float(a[key])
        r.update(case_id=a['case_id'],design_id=a['design_id'],polarization=a['polarization'],
                 bottom_width_um=15.,A_loss_total=float(a['FEM_A_loss_total']),
                 A_port=float(a['FEM_A_port']),A_film=float(a['FEM_A_film']),A_Ti_H=float(a['FEM_A_Ti_H']),
                 ray_A_total=float(a['ray_A_total']),quality_class=a['quality_class'],
                 method_review_failed=a['quality_class']=='method_review_failed_sensitivity')
        rows.append(r)
    assert len(rows)==88 and len({r['case_id'] for r in rows})==88
    baseline_ids={r['case_id'] for r in dataset['baseline_rows']}
    legacy=[r for r in rows if r['case_id'] in baseline_ids]
    pools={'primary_exploratory':[r for r in rows if not r['method_review_failed']],
           'method_checked_only':[r for r in rows if r['quality_class']=='method_checked_exploratory'],
           'legacy_all_warnings_sensitivity':rows}
    # Fixed formulas/hyperparameters, not fitted using held-out outcomes.
    tasks=[('condition_legacy26',legacy,'condition')]+[(k,v,'pooled') for k,v in pools.items()]
    oof=[];metric_rows=[];model_index=[];leak_checks=[]
    for pool,data,mode in tasks:
        conds=sorted({condition(r) for r in data}) if mode=='condition' else ['pooled']
        for c in conds:
            subset=[r for r in data if mode=='pooled' or condition(r)==c]
            groups=sorted({r['design_id'] for r in subset})
            methods=['mean','ridge','gp'] if mode=='condition' else METHODS
            for method in methods:
                method_records=[]
                for group in groups:
                    train=[r for r in subset if r['design_id']!=group]
                    test=[r for r in subset if r['design_id']==group]
                    assert train and not any(r['design_id']==group for r in train)
                    model=fit(train,method,mode)
                    p,s=predict(model,test)
                    assert np.isfinite(p).all()
                    leak_checks.append({'pool':pool,'condition':c,'method':method,'heldout':group,
                                        'train_designs':model['train_designs'],'disjoint':True})
                    for i,r in enumerate(test):
                        item={'pool':pool,'condition':c,'method':method,'case_id':r['case_id'],
                              'heldout_design':group,'train_designs':model['train_designs'],
                              'quality_class':r['quality_class'],'observed':r['A_loss_total'],
                              'predicted':float(p[i]),'std_uncalibrated':float(s[i]) if s is not None else None,
                              'theta_deg':r['theta_deg'],'wavelength_um':r['wavelength_um']}
                        method_records.append(item);oof.append(item)
                stat=metrics(method_records);stat.update(pool=pool,condition=c,method=method)
                metric_rows.append(stat)
                final=fit(subset,method,mode)
                name=f'{pool}__{c}__{method}.joblib';joblib.dump(final,HERE/'models'/name)
                model_index.append({'pool':pool,'condition':c,'method':method,'file':'models/'+name,
                                   'n_rows':len(subset),'n_designs':len(groups),'warnings':final['warnings']})
        print(pool,'complete',round(time.monotonic()-start,2),flush=True)
    dump('dataset.json',{'rows':rows,'pools':{k:[r['case_id'] for r in v] for k,v in pools.items()},
                         'quality_policy':'Method-review failures excluded from primary, retained in separate sensitivity pool. Unreviewed rows remain explicitly exploratory.'})
    dump('heldout_predictions.json',oof);dump('metrics.json',metric_rows)
    dump('model_index.json',model_index);dump('leakage_checks.json',leak_checks)
    dump('manifest.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
         'training_data_as_of':'2026-11-04','prospective_widthfill_as_of':'2026-11-05',
         'python':platform.python_version(),'versions':{'python':platform.python_version(),'numpy':np.__version__,
         'scipy':scipy.__version__,'sklearn':sklearn.__version__,'joblib':joblib.__version__},
         'elapsed_s':time.monotonic()-start,'max_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
         'legacy_rows':88,'independent_training_designs':8,'new_heldout_designs':2,'models_saved':len(model_index),
         'primary_rows':len(pools['primary_exploratory']),'method_checked_rows':len(pools['method_checked_only']),
         'failure_sensitivity_rows':sum(r['method_review_failed'] for r in rows),
         'source_files':[{'path':str(p.relative_to(HERE)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in evidence],
         'validation':'Leave entire H/b/t design out across all wavelengths/angles/polarizations. All scalers fit within training fold.',
         'limitations':['Only eight original designs; b constant at 15um in training.',
           'Primary includes explicitly unreviewed exploratory labels; sensitivity excludes/reincludes method failures separately.',
           'No full 8-13um spectrum/cone target; predictions are point directional absorptivity.',
           'Cross-validation is developmental, not an independent prospective performance guarantee.',
           'GP uncertainty not calibrated; ray model neglects diffraction and whole-cavity coherence.'],
         'server_used':False,'new_FEM_solves':0})

if __name__=='__main__':main()
