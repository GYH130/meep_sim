import ast
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from datetime import datetime, timezone
from unittest.mock import patch

from ti2d.queue import Queue, SolverLedger, aggregate_gate
from ti2d.common import implementation_digest
from ti2d.session_guard import session_live_members, stop_owned_session


class ExpiryTests(unittest.TestCase):
    def test_absolute_deadline_beats_unused_original_budget(self):
        with tempfile.TemporaryDirectory() as d:
            q=object.__new__(Queue)
            q.ledger=SolverLedger(Path(d)/'ledger.json')
            q.budget=dict(total_active_seconds=172800,hard_deadline_epoch=1020)
            with patch('ti2d.queue.time.time',return_value=1000):
                self.assertEqual(q.remaining_seconds(),20)
            with patch('ti2d.queue.time.time',return_value=1030):
                self.assertEqual(q.remaining_seconds(),0)

    def test_partial_suite_never_unlocks_production(self):
        manifest=dict(qualification_subset_only=True,production_gate_task_ids=['x'],
            gate_comparisons=[dict(kind=k,ids=['x'],tolerance=.01) for k in ['mesh','mirror','strict']])
        result=dict(status='QUALIFIED',metrics=dict(R=.9,T=0,A_flux=.1))
        self.assertFalse(aggregate_gate(manifest,{'x':result})['passed'])

    def test_solver_option_and_cache_identity_are_explicit(self):
        root=Path(__file__).resolve().parents[1]
        tree=ast.parse((root/'ti2d/solver.py').read_text())
        simulations=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and
                     isinstance(n.func,ast.Attribute) and n.func.attr=='Simulation']
        self.assertEqual(len(simulations),1)
        self.assertTrue(any(k.arg=='eps_averaging' for k in simulations[0].keywords))
        self.assertIn("eps_averaging=c['eps_averaging']",(root/'ti2d/solver.py').read_text())
        self.assertNotEqual(implementation_digest(),'24d09c50dcf6c0cb9271947dbcb04f03b62d9d0500bd7bace34e6909d54834c1')

    def test_guard_completed_session(self):
        class Child:
            pid=123
            returncode=0
            def poll(self): return 0
        with patch('ti2d.session_guard.session_live_members',return_value=[]):
            self.assertTrue(stop_owned_session(Child(),0,lambda:None)['stopped'])

    def test_guard_does_not_claim_orphan_rank_stopped(self):
        class Child:
            pid=123
            returncode=0
            def poll(self): return 0
        with patch('ti2d.session_guard.session_live_members',return_value=[{'pid':124}]):
            self.assertFalse(stop_owned_session(Child(),0,lambda:None)['stopped'])

    def test_creator_has_no_meep_import(self):
        # Keep the frozen H creator regression, not the unrelated new K creator.
        root=Path(__file__).resolve().parents[2] / 'expiry_campaign'
        spec=importlib.util.spec_from_file_location('expiry_create',root/'create_campaign.py')
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(module.DEADLINE,datetime(2026,9,23,8,52,44,tzinfo=timezone.utc).timestamp())
        self.assertEqual(module.RUN.name,'expiry_20260922_H')


if __name__=='__main__': unittest.main()
