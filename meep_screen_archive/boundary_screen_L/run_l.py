"""Finite homogeneous-boundary screen. Never unlocks production or changes K."""
from pathlib import Path
import json
import math
import os
import signal
import sys
import time

from ti2d.common import atomic_json, digest_file, utcnow
from ti2d.queue import Queue, SolverLedger, FileLock
from prepare_l import ROOT, RUN, CONTROL_ID, CANDIDATE_IDS


def control_audit(result, proof):
    reference = result.get('reference', {})
    pseudo = reference.get('pseudo_reflection')
    power = reference.get('input_power')
    width = result.get('case_config', {}).get('geometry', {}).get('period_um')
    expected = proof['pseudo_reflection']
    expected_power = proof['input_power_per_width']
    valid = all(isinstance(x, (int,float)) and math.isfinite(x) for x in
                (pseudo,power,width,expected,expected_power)) and width > 0 and expected_power > 0
    checks = reference.get('checks', {})
    physical = all(checks.get(k) is True for k in
                   ('stop','downward_flux','direction','plane_closure','empty_transmission'))
    passed = (valid and physical and result.get('status') == 'NUMERICAL_FAILED'
              and result.get('failure_stage') == 'reference'
              and checks.get('empty_pseudo_reflection') is False
              and abs(pseudo-expected) <= 1e-5
              and abs((power/width)/expected_power-1) <= 1e-3)
    return dict(reproduced=bool(passed), pseudo_reflection=pseudo,
                expected_pseudo_reflection=expected, absolute_tolerance=1e-5,
                input_power_per_width=power/width if valid else None,
                expected_power_per_width=expected_power, relative_power_tolerance=1e-3,
                qualification_unchanged=True,
                note='Reproducing a known failure validates this diagnostic comparison, not the boundary.')


def candidate_pass(result):
    pseudo = result.get('metrics', {}).get('reference_pseudo_reflection')
    return (result.get('status') == 'QUALIFIED' and isinstance(pseudo,(int,float))
            and math.isfinite(pseudo) and pseudo < 5e-4)


class LQueue(Queue):
    def __init__(self):
        super().__init__(RUN/'manifest.json', ('validation',))
        for field, base in (('scripts_sha256', ROOT), ('inputs_sha256', RUN),
                            ('future_configs_sha256', RUN), ('material_files_sha256', ROOT)):
            records = self.manifest.get(field)
            if not records:
                raise RuntimeError('MISSING_FROZEN_INPUTS:' + field)
            for name, expected in records.items():
                path = base / name
                if (path.is_symlink() or not path.resolve().is_relative_to(base.resolve())
                        or not path.is_file() or digest_file(path) != expected):
                    raise RuntimeError('FROZEN_INPUT_CHANGED:' + name)
        self.audit = None
        self.selected = None

    def execute(self, task, attempt, admission=None):
        original = self.budget['single_solver_seconds']
        self.budget['single_solver_seconds'] = min(1200., original)
        self.statuses[task['id']] = {'status':'RUNNING'}
        try:
            return super().execute(task, attempt, admission)
        finally:
            self.budget['single_solver_seconds'] = original

    def checkpoint(self, final=False):
        from export_l import main as export
        folder = RUN/'reports'
        folder.mkdir(exist_ok=True)
        rows=[]
        for task in self.manifest['tasks']:
            result = self.results.get(task['id'], {})
            reference = result.get('reference', {})
            metrics = result.get('metrics', {})
            rows.append(dict(id=task['id'], status=self.statuses.get(task['id'],{}).get('status','NOT_RUN'),
                failure_stage=result.get('failure_stage'),
                reference_pseudo_reflection=metrics.get('reference_pseudo_reflection',reference.get('pseudo_reflection')),
                R=metrics.get('R'), T=metrics.get('T'), A_flux=metrics.get('A_flux'), A_vol=metrics.get('A_vol'),
                elapsed_s=result.get('queue_process',{}).get('elapsed_s')))
        summary = dict(updated_utc=utcnow(),final=final,control_audit=self.audit,selected_candidate=self.selected,
            conditions=rows,usage=self.ledger.as_dict(),production_qualified=False,
            full_width_validation='DEFERRED_RENTAL_WINDOW',
            scope='Uniform air/slab numerical width 2 um; not a machining-period study or slotted-geometry qualification.',
            limits='Selected profile remains unverified for the full 50 um cell and grazing slot diffraction orders.',
            ai_usage={'available':False,'rmb':None,'queue_ai_calls':0})
        atomic_json(folder/'screen_summary.json',summary)
        text=['# L — homogeneous boundary screening only','',summary['scope'],'',summary['limits'],'',
              'Full 50 um reference/flat/slot: DEFERRED_RENTAL_WINDOW. No production gate unlocked.','',
              '| Condition | Solver status | Background reflected power ratio |', '|---|---|---|']
        for row in rows:
            text.append('| '+row['id']+' | '+row['status']+' | '+str(row['reference_pseudo_reflection'])+' |')
        text += ['', 'Selected candidate: '+str(self.selected),
                 'Original pseudo-reflection acceptance remains <0.001; screening additionally requires <0.0005.',
                 'Control must reproduce K at absolute difference <=0.00001 and per-width incident power difference <=0.1%.',
                 'Maximum 3 conditions, 20 minutes each, total solver union 1 hour; no retries.',
                 '2026-09-24 20:30 China solver cutoff; approximately21:53 rental expiry.',
                 'Raw arrays are backed up separately; only sanitized compact evidence goes to the archive branch.']
        (folder/'L_HANDOFF.md').write_text('\n'.join(text)+'\n')
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            selected=[r for r in rows if isinstance(r['reference_pseudo_reflection'],(int,float))]
            if selected:
                fig,ax=plt.subplots(figsize=(8,4),constrained_layout=True)
                ax.bar(range(len(selected)),[r['reference_pseudo_reflection'] for r in selected])
                ax.set_xticks(range(len(selected)),[r['id'].replace('L_narrow_','') for r in selected])
                ax.axhline(.001,color='red',ls='--',label='Original gate')
                ax.axhline(.0005,color='green',ls=':',label='Screening margin')
                ax.set(ylabel='Backward / forward power',title='Uniform planar boundary screen — not slot emissivity')
                ax.legend(); fig.savefig(folder/'boundary_screen.png',dpi=150); plt.close(fig)
        except Exception as exc:
            atomic_json(folder/'figure_error.json',{'error':str(exc)})
        return export(final=final)

    def run_screen(self):
        with FileLock(RUN/'queue.lock'), FileLock(self.ledger.path.with_suffix('.lock')):
            self.ledger=SolverLedger(self.ledger.path)
            if any(r.get('end') is None for r in self.ledger.data['intervals']):
                raise RuntimeError('OPEN_PRIOR_INTERVAL')
            if (RUN/'lock_manifest.json').exists():
                raise RuntimeError('SINGLE_USE_SCREEN_NO_AUTOMATIC_RETRY')
            atomic_json(RUN/'lock_manifest.json',self.locked)
            signal.signal(signal.SIGTERM,self.cancel)
            signal.signal(signal.SIGINT,self.cancel)
            self.statuses={t['id']:{'status':'NOT_RUN'} for t in self.manifest['tasks']}
            proof=json.loads((RUN/'inputs/K_reference_scalars.json').read_text())
            state='SCREEN_INCOMPLETE'
            try:
                for task in self.manifest['tasks']:
                    if self.cancelled or self.remaining_seconds()<180:
                        state='SCREEN_BUDGET_OR_INTERRUPT'; break
                    if task['id'] != CONTROL_ID and not (self.audit or {}).get('reproduced'):
                        state='CONTROL_NOT_REPRODUCED'; break
                    result, directory=self.execute(task,1)
                    self.results[task['id']]=result
                    self.statuses[task['id']]={'status':result['status'],'result_path':str(directory/'result.json')}
                    if task['id']==CONTROL_ID:
                        self.audit=control_audit(result,proof)
                        if not self.audit['reproduced']: state='CONTROL_NOT_REPRODUCED'
                    elif candidate_pass(result):
                        self.selected=task['id']; state='SCREEN_CANDIDATE_FOUND_FULL_VALIDATION_NOT_RUN'
                    elif result['status'] not in ('NUMERICAL_FAILED','QUALIFIED'):
                        state='SCREEN_ABORTED_NO_RETRY'
                    else:
                        state='SCREEN_NO_CANDIDATE_PASSED'
                    self.snapshot('CHECKPOINT')
                    self.checkpoint()
                    if self.selected or state in ('CONTROL_NOT_REPRODUCED','SCREEN_ABORTED_NO_RETRY'):
                        break
                self.snapshot(state)
            except Exception as exc:
                self.snapshot('QUEUE_ERROR',message=str(exc))
                raise
            finally:
                self.checkpoint(final=True)


if __name__=='__main__':
    if sys.platform!='linux' or str(ROOT)!='/root/autodl-tmp/meep_sim/meep_screen/boundary_screen_L':
        raise RuntimeError('SERVER_ONLY_NO_LOCAL_FDTD')
    os.chdir(ROOT)
    for key,value in {'MPLCONFIGDIR':ROOT/'cache/matplotlib','XDG_CACHE_HOME':ROOT/'cache','TMPDIR':ROOT/'tmp'}.items():
        value.mkdir(parents=True,exist_ok=True); os.environ[key]=str(value)
    LQueue().run_screen()
