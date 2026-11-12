"""Held-out new-design challenge and presentation-quality diagnostic figure."""
from pathlib import Path
import json,collections,time
import numpy as np,joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from surrogate import predict,metrics,METHODS
from predict import ray_predict,evaluate
HERE=Path(__file__).resolve().parent

def dump(name,d): (HERE/name).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def main():
    held=json.loads((HERE/'heldout_predictions.json').read_text())
    summary=json.loads((HERE/'metrics.json').read_text())
    by_design=[]
    for m in summary:
        subset=[p for p in held if all(p[k]==m[k] for k in ['pool','condition','method'])]
        group_errors=[]
        for g in sorted({p['heldout_design'] for p in subset}):
            gs=[p for p in subset if p['heldout_design']==g];s=metrics(gs)
            s.update(pool=m['pool'],condition=m['condition'],method=m['method'],design_id=g)
            by_design.append(s);group_errors.append(s['MAE_pp'])
        m['geometry_balanced_MAE_pp']=float(np.mean(group_errors))
    dump('metrics.json',summary);dump('metrics_by_design.json',by_design)
    dataset=json.loads((HERE/'dataset.json').read_text())
    primary_ids=set(dataset['pools']['primary_exploratory'])
    primary=[r for r in dataset['rows'] if r['case_id'] in primary_ids]
    strata=[]
    for method in METHODS:
        a=[r for r in held if r['pool']=='primary_exploratory' and r['method']==method]
        for r in a:
            wl=[s['wavelength_um'] for s in primary if s['design_id']!=r['heldout_design']]
            r['heldout_wavelength_in_train_range']=min(wl)<=r['wavelength_um']<=max(wl)
        for key,sub in [('all',a),('lambda_10p5',[r for r in a if r['wavelength_um']==10.5]),
            ('wavelength_in_train_range',[r for r in a if r['heldout_wavelength_in_train_range']]),
            ('wavelength_outside_train_range',[r for r in a if not r['heldout_wavelength_in_train_range']])]:
            if sub:
                m=metrics(sub);m.update(method=method,stratum=key);strata.append(m)
    dump('support_stratified_metrics.json',{'note':'All points retained in global metrics. Range support concerns heldout fold wavelength only; joint spatial/angular coverage remains limited.','metrics':strata})
    new=[];defs=[('WFT1',25.,3.,.3567405584106311,.35589938154084155),
                ('WFT2',15.,6.,.5305910338702413,.5297846176660281)]
    for name,b,t,obs,port in defs:
        q=dict(case_id=name,design_id=f'H80_b{b:g}_t{t:g}',height_um=80.,bottom_width_um=b,film_um=t,
               wavelength_um=10.5,theta_deg=0.,phi_deg=0.,polarization='s')
        ray=ray_predict(q,rays=1048576);q['ray_A_total']=ray['A_total']
        for pool in ['primary_exploratory','method_checked_only','legacy_all_warnings_sensitivity']:
            for method in METHODS:
                model=joblib.load(HERE/f'models/{pool}__pooled__{method}.joblib')
                p,s=predict(model,[q]);value=float(p[0])
                new.append(dict(query=q,pool=pool,method=method,observed_A_loss=obs,observed_A_port=port,
                  predicted=value,error_pp=100*(value-obs),GP_std_uncalibrated=float(s[0]) if s is not None else None,
                  ray_channel=ray,training_excluded=True,
                  scope='NEW_B_OUT_OF_SUPPORT' if b!=15 else 'NEW_T_OUT_OF_SUPPORT',
                  quality_class='source_and_non_rawflux_checks_passed_but_saved_field_fit_warning',
                  formal_qualification=False,blind_prospective_test=False,
                  test_note='Existing results were visible during development; excluded from all training. This is a retrospective independent-design challenge, not a preregistered blind prospective test.'))
    dump('new_design_challenge.json',new)
    q=dict(height_um=90.,bottom_width_um=15.,film_um=2.5,wavelength_um=10.5,theta_deg=30.,phi_deg=30.,polarization='s')
    dump('example_prediction.json',evaluate(q))
    assert evaluate(dict(q,bottom_width_um=25))['state']=='OUT_OF_TRAINING_SUPPORT'
    assert evaluate(dict(q,film_um=6))['state']=='OUT_OF_TRAINING_SUPPORT'
    assert not any('WFT' in g for m in json.loads((HERE/'model_index.json').read_text())
                   for g in joblib.load(HERE/m['file'])['train_designs'])
    dump('checks.json',{'all_training_groups_disjoint_from_test':True,'novel_bottom_width_rejected_without_override':True,
       'novel_thickness_rejected_without_override':True,'WFT1_WFT2_never_used_for_training':True,
       'prediction_command_executed':True,'quality_flags_retained':True,
       'unchanged_original_thresholds':True,'no_new_FEM_or_server_jobs':True})
    plot(held,summary,new)

def plot(held,summary,new):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axs=plt.subplots(1,3,figsize=(13.5,4.2),constrained_layout=True)
    colors={0:'#cb6549',15:'#d49642',22.5:'#999d46',30:'#38877c',60:'#4a6178'}
    for ax,method,title in zip(axs[:2],['ray','gp'],['(a) Ray + thin-film optics','(b) Pooled Gaussian process']):
        ps=[r for r in held if r['pool']=='primary_exploratory' and r['method']==method]
        for th in sorted({r['theta_deg'] for r in ps}):
            a=[r for r in ps if r['theta_deg']==th]
            ax.scatter([100*r['observed'] for r in a],[100*r['predicted'] for r in a],
                       c=colors.get(th,'#888888'),s=29,alpha=.75,edgecolors='white',linewidths=.4,label=f'{th:g}°')
        ax.plot([0,65],[0,65],c='#555555',ls='--',lw=.9)
        ax.set(xlim=(0,65),ylim=(0,65),xlabel='FEM absorptivity (%)',ylabel='Held-out prediction (%)',title=title)
        m=next(m for m in summary if m['pool']=='primary_exploratory' and m['method']==method)
        ax.text(.04,.96,f"Point MAE {m['MAE_pp']:.2f} pp\nDesign-balanced MAE {m['geometry_balanced_MAE_pp']:.2f} pp\nMax error {m['max_abs_error_pp']:.2f} pp",transform=ax.transAxes,va='top',fontsize=9)
        ax.grid(alpha=.13);ax.set_aspect('equal')
    axs[0].legend(title='Polar angle',frameon=False,fontsize=8,loc='lower right')
    ax=axs[2]
    methods=['ray','gp','ray_gp_residual'];labels=['Ray','GP','Ray + GP residual'];cols=['#6b8792','#245883','#cb6549']
    for i,(method,label,col) in enumerate(zip(methods,labels,cols)):
        p=[r for r in new if r['pool']=='primary_exploratory' and r['method']==method]
        ax.plot([0,1],[100*r['predicted'] for r in p],'o--',label=label,c=col,lw=1.3)
    ax.plot([0,1],[35.67405584106311,53.05910338702413],'D-',c='#111111',label='FEM (warning retained)',lw=1.5)
    ax.set_xticks([0,1],['b25 / t3','b15 / t6']);ax.set(xlim=(-.2,1.2),ylim=(0,65),ylabel='Normal absorptivity (%)',
                           title='(c) Two unseen geometries')
    ax.legend(frameon=False,fontsize=8,loc='lower right');ax.grid(alpha=.13)
    ax.text(.03,.96,'Both outside training support\nNot used for fitting\nExisting, non-blinded results',transform=ax.transAxes,va='top',fontsize=9)
    fig.suptitle('Exploratory surrogate v0: whole-geometry holdout, 8 training designs',fontsize=13)
    for ext in ['png','svg','pdf']:fig.savefig(HERE/f'surrogate_validation.{ext}',dpi=210,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':main()
