import hashlib
import io
import json
import unittest

import numpy as np

from recover_tail_analysis import boundary_comparison, checked_npz, load_waveforms, regular_selection, resolve_status
from tail_analysis import replay_audit, spectrum


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.t=np.arange(0,640,.25)
        self.fields=np.column_stack([np.exp(-2j*np.pi*.1*self.t),np.zeros(len(self.t))])

    def test_uniform_grid_and_complex_signed_frequency(self):
        keep,audit=regular_selection(self.t,self.fields,.25)
        self.assertEqual(len(keep),len(self.t))
        self.assertEqual(audit['excluded_samples'],[])
        _,_,detail=spectrum(self.t[keep],self.fields[keep])
        self.assertAlmostEqual(detail['peaks'][0]['signed_frequency'],-.1)

    def test_terminal_partial_is_explicit_and_original_untouched(self):
        t=np.r_[self.t,self.t[-1]+.125]
        f=np.vstack([self.fields,self.fields[-1]])
        oldt=t.copy();oldf=f.copy()
        keep,audit=regular_selection(t,f,.25)
        self.assertEqual(len(keep),len(t)-1)
        self.assertEqual(audit['excluded_samples'][0]['actual_dt'],.125)
        self.assertFalse(audit['resampled'])
        np.testing.assert_array_equal(t,oldt);np.testing.assert_array_equal(f,oldf)

    def test_interior_gap_rejected_even_with_terminal_partial(self):
        t=np.delete(self.t,200);f=np.delete(self.fields,200,axis=0)
        with self.assertRaisesRegex(ValueError,'TIME_GRID_DEFECT'):
            regular_selection(t,f,.25)
        with self.assertRaisesRegex(ValueError,'TIME_GRID_DEFECT'):
            regular_selection(np.r_[t,t[-1]+.125],np.vstack([f,f[-1]]),.25)

    def test_duplicate_and_reverse_rejected(self):
        t=self.t.copy();t[200]=t[199]
        with self.assertRaisesRegex(ValueError,'DUPLICATE_OR_REVERSED'):
            regular_selection(t,self.fields,.25)
        t[200]=t[198]
        with self.assertRaisesRegex(ValueError,'DUPLICATE_OR_REVERSED'):
            regular_selection(t,self.fields,.25)

    def test_nonfinite_times_fields_and_terminal_sample_rejected(self):
        for index in (20,-1):
            t=self.t.copy();t[index]=np.nan
            with self.assertRaisesRegex(ValueError,'NONFINITE'):
                regular_selection(t,self.fields,.25)
            f=self.fields.copy();f[index,0]=np.inf
            with self.assertRaisesRegex(ValueError,'NONFINITE'):
                regular_selection(self.t,f,.25)

    def test_checksum_and_array_finiteness_rejected(self):
        buf=io.BytesIO();np.savez_compressed(buf,t=self.t,fields=self.fields)
        payload=buf.getvalue();sha=hashlib.sha256(payload).hexdigest()
        self.assertIn('t',checked_npz(payload,sha))
        with self.assertRaisesRegex(ValueError,'CHECKSUM'):
            checked_npz(payload,'0'*64)
        buf=io.BytesIO();np.savez_compressed(buf,t=np.array([np.nan]))
        payload=buf.getvalue()
        with self.assertRaisesRegex(ValueError,'NONFINITE'):
            checked_npz(payload,hashlib.sha256(payload).hexdigest())

    def test_final_full_missing_step_cannot_be_excluded(self):
        t=self.t.copy();t[-1]+=.25
        with self.assertRaisesRegex(ValueError,'TIME_GRID_DEFECT'):
            regular_selection(t,self.fields,.25)

    def test_replay_gate_not_relaxed_by_report_recovery(self):
        h=[dict(time=float(t),probe_power={'back_air':.001}) for t in range(2,1201,2)]
        self.assertFalse(replay_audit(h,h[:457],1.)['passed'])
        self.assertTrue(replay_audit(h,h,1.)['passed'])

    def test_separate_terminal_artifact_must_match_result_and_hash(self):
        def pack(**kw):
            b=io.BytesIO();np.savez_compressed(b,**kw);return b.getvalue()
        t=np.array([.25,.5,.75]);v=np.zeros((3,1,3),complex)
        waveform=pack(point_t=t,point_fields=v,line_t=np.array([]))
        terminal=pack(point_t=np.array([.875]),point_fields=np.zeros((1,1,3),complex))
        files={'arm/waveforms_0000.npz':waveform,'arm/terminal_sample.npz':terminal,
               'arm/terminal_sample.json':json.dumps(dict(file='terminal_sample.npz',
                    sha256=hashlib.sha256(terminal).hexdigest(),time=.875)).encode()}
        manifest=[dict(file='waveforms_0000.npz',sha256=hashlib.sha256(waveform).hexdigest(),first_time=.25,last_time=.75)]
        metadata=dict(probe_names=['front_air'],point_component_names=['Ex','Ey','Hz'],point_dt=.25)
        _,_,_,audit=load_waveforms(files.__getitem__,'arm',manifest,metadata,dict(simulation_time=.875))
        self.assertTrue(audit['separate_terminal_sample']['excluded_from_fft'])
        with self.assertRaisesRegex(ValueError,'TERMINAL_SAMPLE_RESULT_MISMATCH'):
            load_waveforms(files.__getitem__,'arm',manifest,metadata,dict(simulation_time=.9))
        files['arm/terminal_sample.npz']=terminal+b'corrupt'
        with self.assertRaisesRegex(ValueError,'CHECKSUM'):
            load_waveforms(files.__getitem__,'arm',manifest,metadata,dict(simulation_time=.875))

    def comparison_inputs(self):
        t=np.arange(.25,1200.25,.25)
        values=np.zeros((len(t),2,3),complex)
        values[:,:,0]=np.exp(-2j*np.pi*.1*t)[:,None]
        arms={name:dict(status='DIAGNOSTIC_COMPLETE',worker_status='DIAGNOSTIC_COMPLETE') for name in ('pml6_control','pml12_test')}
        waveforms={name:dict(point_t=t.copy(),point_fields=values.copy()*scale,
                             point_dt=.25,probe_names=['front_air','back_air'],source_end=525.)
                   for name,scale in (('pml6_control',1.),('pml12_test',.5))}
        return arms,waveforms

    def test_matched_boundary_window_stats_and_signed_spectrum(self):
        arms,waveforms=self.comparison_inputs()
        result=boundary_comparison(arms,waveforms,True,2.)
        self.assertEqual(result['interval'],[1000.,1200.])
        self.assertEqual(result['samples'],801)
        self.assertAlmostEqual(result['control_mean_intensity_ratio'],.5)
        self.assertAlmostEqual(result['thick_over_control'],.25)
        self.assertAlmostEqual(result['thick_over_control_max'],.25)
        self.assertFalse(result['qualification'])
        self.assertAlmostEqual(result['matched_window_spectra']['pml6_control']['back_air']['peaks'][0]['signed_frequency'],-.1,places=3)

    def test_boundary_mismatched_grids_rejected(self):
        arms,waveforms=self.comparison_inputs()
        waveforms['pml12_test']['point_t']+=.125
        with self.assertRaisesRegex(ValueError,'COMPARISON_TIME_GRIDS_MISMATCH'):
            boundary_comparison(arms,waveforms,True,1.)

    def test_boundary_not_evaluated_without_replay_or_completed_arms(self):
        arms,waveforms=self.comparison_inputs()
        self.assertIsNone(boundary_comparison(arms,waveforms,False,1.))
        arms['pml12_test']['status']='TIME_LIMIT_UNQUALIFIED'
        self.assertIsNone(boundary_comparison(arms,waveforms,True,1.))

    def test_supervisor_failure_overrides_worker_termination_label(self):
        for cause in ('MEMORY_LIMIT','INTERRUPTED','CRASHED'):
            result=resolve_status({'status':'TIME_LIMIT_UNQUALIFIED'},{'status':cause})
            self.assertEqual(result['status'],cause)
            self.assertEqual(result['worker_status'],'TIME_LIMIT_UNQUALIFIED')
            self.assertEqual(result['status_authority'],'process.json')

    def test_process_status_without_worker_is_not_not_run(self):
        result=resolve_status(None,{'status':'MEMORY_LIMIT'})
        self.assertEqual(result['status'],'MEMORY_LIMIT')
        self.assertIsNone(result['worker_status'])
        self.assertEqual(resolve_status()['status'],'NOT_RUN_OR_NO_RESULT')

    def test_boundary_requires_worker_and_supervisor_complete(self):
        arms,waveforms=self.comparison_inputs()
        arms['pml12_test']=resolve_status({'status':'TIME_LIMIT_UNQUALIFIED'},{'status':'DIAGNOSTIC_COMPLETE'})
        self.assertIsNone(boundary_comparison(arms,waveforms,True,1.))
        arms['pml12_test']=resolve_status({'status':'DIAGNOSTIC_COMPLETE'},{'status':'INTERRUPTED'})
        self.assertIsNone(boundary_comparison(arms,waveforms,True,1.))


if __name__=='__main__':unittest.main()
