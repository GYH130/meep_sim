import unittest
from unittest.mock import Mock
import numpy as np
from tail_analysis import replay_audit, spectrum
from run_tail_diagnostic import allowance, DEADLINE
from tail_probe import probe_positions
from ti2d.solver import layout


class TailTests(unittest.TestCase):
    def test_complex_tone_frequency_and_sampling(self):
        t=np.arange(0,640,.25)
        e=np.exp(-2j*np.pi*.1*t)
        _,_,s=spectrum(t,np.column_stack([e,np.zeros_like(e)]))
        self.assertAlmostEqual(s['peaks'][0]['signed_frequency'],-.1)
        self.assertEqual(s['nyquist'],2.)
        with self.assertRaises(ValueError): spectrum(np.r_[t[:100],t[101:]],np.column_stack([e[:-1],e[:-1]]))

    def test_control_requires_matching_post_source_evidence(self):
        h=[dict(time=float(t),probe_power={'back_air':.001}) for t in range(2,1201,2)]
        self.assertTrue(replay_audit(h,h,1.)['passed'])
        self.assertFalse(replay_audit(h,h[:200],1.)['passed'])
        wrong=[dict(r,probe_power={'back_air':.01}) for r in h]
        self.assertFalse(replay_audit(h,wrong,1.)['passed'])

    def test_budget_obeys_all_caps(self):
        ledger=Mock(); ledger.total_active_seconds.return_value=64000
        plan={'ledger_active_before':63000}
        self.assertEqual(allowance(plan,ledger,DEADLINE-500),380)
        self.assertEqual(allowance(plan,ledger,DEADLINE+1),0)
        ledger.total_active_seconds.return_value=80950
        self.assertEqual(allowance(plan,ledger,DEADLINE-20000),50)
        ledger.total_active_seconds.return_value=172810
        self.assertEqual(allowance(plan,ledger,DEADLINE-20000),0)

    def test_pml_change_preserves_all_physical_probe_positions(self):
        c=dict(geometry=dict(period_um=50.,axis_length_um=40.,tilt_deg=30.,thickness_um=60.),
               pml_um=6.,air_above_um=12.,air_below_um=12.)
        d=dict(c,pml_um=12.)
        self.assertEqual(probe_positions(c),probe_positions(d))
        self.assertEqual(layout(d)['height']-layout(c)['height'],12.)


if __name__=='__main__': unittest.main()
