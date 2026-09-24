import unittest
from copy import deepcopy
from run_l import control_audit,candidate_pass

class ControlTests(unittest.TestCase):
    def test_expected_failure_is_diagnostic_only(self):
        proof={'pseudo_reflection':.0017166551733512714,'input_power_per_width':732.155420649526}
        r={'status':'NUMERICAL_FAILED','failure_stage':'reference','case_config':{'geometry':{'period_um':2}},
           'reference':{'pseudo_reflection':proof['pseudo_reflection'],'input_power':2*proof['input_power_per_width'],
                        'checks':{k:True for k in ('stop','downward_flux','direction','plane_closure','empty_transmission')}}}
        r['reference']['checks']['empty_pseudo_reflection']=False
        self.assertTrue(control_audit(r,proof)['reproduced'])
        changed=deepcopy(r);changed['reference']['pseudo_reflection']+=.00002
        self.assertFalse(control_audit(changed,proof)['reproduced'])
        changed=deepcopy(r);changed['reference']['input_power']*=1.01
        self.assertFalse(control_audit(changed,proof)['reproduced'])
        self.assertFalse(candidate_pass(r))

    def test_margin_does_not_replace_other_gates(self):
        self.assertFalse(candidate_pass({'status':'NUMERICAL_FAILED','metrics':{'reference_pseudo_reflection':1e-5}}))
        self.assertFalse(candidate_pass({'status':'QUALIFIED','metrics':{'reference_pseudo_reflection':.0007}}))
        self.assertTrue(candidate_pass({'status':'QUALIFIED','metrics':{'reference_pseudo_reflection':.0001}}))

if __name__=='__main__': unittest.main()
