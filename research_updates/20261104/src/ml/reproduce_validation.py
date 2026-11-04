"""Local-only, portable leave-one-whole-design-out surrogate validation.

Uses preserved scientific fitting functions. Does not open MPH, contact servers,
load pickled models, run a solver, or claim formal physical qualification.
"""
from pathlib import Path
import json
import sys
import numpy as np
import model_helpers as helper

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'reproduced'
METHODS = ['mean', 'linear_ridge', 'gp_matern']

def run():
    dataset = json.loads((ROOT/'data/ml/dataset.json').read_text())
    assert (len(dataset['records_all']), len(dataset['rows']),len(dataset['design_ids'])) == (89,86,8)
    assert len({(r['design_id'],helper.cond(r)) for r in dataset['baseline_rows']}) == 26
    outputs=[]
    pools={'legacy_exploratory': [r for r in dataset['baseline_rows'] if r['exploratory_fit_eligible']],
           'method_checked_sensitivity': [r for r in dataset['baseline_rows'] if r['method_checked_fit_eligible']]}
    for pool, rows in pools.items():
        predictions=[]; metrics=[]; warning_rows=[]
        for condition in sorted({helper.cond(r) for r in rows}):
            subset=[r for r in rows if helper.cond(r)==condition]
            designs=sorted({r['design_id'] for r in subset})
            if len(designs)<3:
                continue
            for label in helper.LABELS:
                for method in METHODS:
                    part=[]
                    for heldout in designs:
                        train=[r for r in rows if r['design_id']!=heldout and helper.cond(r)==condition]
                        test=[r for r in subset if r['design_id']==heldout]
                        assert len(test)==1 and heldout not in {r['design_id'] for r in train}
                        model=helper.fit(train,label,method)
                        mean,std=helper.predict(model,helper.features(test))
                        item={'pool':pool,'condition':condition,'label':label,'method':method,
                              'design_id':heldout,'train_designs':model['train_designs'],
                              'observed':test[0][label],'predicted':float(mean[0]),
                              'std':float(std[0]) if std is not None else None}
                        part.append(item); predictions.append(item)
                        warning_rows += [{'condition':condition,'label':label,'method':method,
                                          'heldout':heldout,'warnings':model['warnings']}] if model['warnings'] else []
                    error=np.array([r['predicted']-r['observed'] for r in part])
                    metrics.append({'pool':pool,'condition':condition,'label':label,'method':method,
                        'n_conditions':len(part),'n_independent_designs':len(designs),
                        'MAE_pp':float(np.mean(abs(error))*100),
                        'RMSE_pp':float(np.sqrt(np.mean(error**2))*100),
                        'max_abs_error_pp':float(max(abs(error))*100),
                        'nominal95_interval_coverage':float(np.mean([abs(p['predicted']-p['observed'])<=1.96*p['std'] for p in part])) if method=='gp_matern' else None})
        expected=json.loads((ROOT/f'data/ml/{pool}_metrics.json').read_text())
        assert len(expected)==len(metrics)
        key=lambda r:(r['condition'],r['label'],r['method'])
        lookup={key(r):r for r in expected}
        max_difference=0.
        for r in metrics:
            e=lookup[key(r)]
            for field in ['MAE_pp','RMSE_pp','max_abs_error_pp']:
                difference=abs(r[field]-e[field]); max_difference=max(max_difference,difference)
                # Optimizers/BLAS versions may change GP last digits; report rather than hide it.
                if difference>0.01:
                    raise AssertionError((pool,key(r),field,r[field],e[field]))
        OUT.mkdir(exist_ok=True)
        for name,value in [('metrics',metrics),('predictions',predictions),('warnings',warning_rows)]:
            (OUT/f'{pool}_{name}.json').write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
        outputs.append({'pool':pool,'models_compared':len(metrics),
                        'heldout_predictions':len(predictions),
                        'max_metric_difference_pp':max_difference,'matches_archive_within_0p01pp':True})
    summary={'state':'LOCAL_VALIDATION_REPRODUCED_NOT_NEW_PHYSICAL_SOLVES',
             'data_rows':86,'independent_designs':8,'formal_qualification':False,
             'checks':outputs,'note':'GP intervals are conditional and uncalibrated; no material/mesh systematic error bound.'}
    (OUT/'validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    run()
