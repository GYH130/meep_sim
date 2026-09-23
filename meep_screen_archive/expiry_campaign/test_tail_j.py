import ast
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch

import run_tail_j as queue
from tail_probe_v2 import probe_positions


class JTests(unittest.TestCase):
    def test_budget_caps_and_removed_calendar_cutoff(self):
        ledger=Mock(); ledger.total_active_seconds.return_value=72000.
        plan={'ledger_active_before':72000.}
        self.assertEqual(queue.allowance(plan,ledger,1790160000.,14400.),14400.)
        self.assertEqual(queue.allowance(plan,ledger,1790160000.,18000.),18000.)
        ledger.total_active_seconds.return_value=104000.
        self.assertEqual(queue.allowance(plan,ledger,0,18000.),400.)
        plan={'ledger_active_before':170000.}
        ledger.total_active_seconds.return_value=172500.
        self.assertEqual(queue.allowance(plan,ledger,0,18000.),300.)
        ledger.total_active_seconds.return_value=173000.
        self.assertEqual(queue.allowance(plan,ledger,0,18000.),0.)

    def test_packaging_preserves_prior_I_bundle(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); run=root/'runs/J'; run.mkdir(parents=True)
            (run/'status.json').write_text('{}')
            exports=root/'exports'; exports.mkdir()
            prior=exports/'I_raw_evidence.tar.gz'; prior.write_bytes(b'old-evidence')
            with patch.object(queue,'ROOT',root), patch.object(queue,'RUN',run), patch.object(queue,'SCRIPTS',()):
                queue.package()
            self.assertEqual(prior.read_bytes(),b'old-evidence')
            self.assertTrue((exports/'J_raw_evidence.tar.gz').is_file())
            self.assertTrue((exports/'J_export.json').is_file())

    def test_probe_physics_unchanged(self):
        root=Path(__file__).parent
        old=ast.parse((root/'tail_probe.py').read_text())
        new=ast.parse((root/'tail_probe_v2.py').read_text())
        def setup(tree):
            main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
            before=[]
            for n in main.body:
                if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='wall' for t in n.targets):break
                before.append(ast.dump(n,include_attributes=False))
            return before
        self.assertEqual(setup(old),setup(new))
        def nested(tree,name):
            return next(ast.dump(n,include_attributes=False) for n in ast.walk(tree)
                        if isinstance(n,ast.FunctionDef) and n.name==name)
        self.assertEqual(nested(old,'collect'),nested(new,'collect'))
        c=dict(geometry=dict(period_um=50.,axis_length_um=40.,tilt_deg=30.,thickness_um=60.),
               pml_um=6.,air_above_um=12.,air_below_um=12.)
        self.assertEqual(probe_positions(c),probe_positions(dict(c,pml_um=12.)))


if __name__=='__main__':unittest.main()
