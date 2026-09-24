"""Queue tests use temporary synthetic evidence, never Meep/FDTD."""
import json
from pathlib import Path
import tempfile
import threading
import time
import unittest
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

from ti2d.queue import (FileLock,SolverLedger,aggregate_gate,atomic_json,
                        case_fingerprint,cgroup_resources,process_tree_rss,
                        reusable_result,seal_result,validate_result,verify_manifest,
                        ledger_path_for_manifest,pair_admission,concurrency_assessment,Queue)
from ti2d.report import paired_metrics,sensitivities,collect


def qualified(value=.25):
    return {"status":"QUALIFIED","stop":{"converged":True,"stop_reason":"CONVERGED"},
            "case_config":{"case":"structure"},
            "metrics":{"R":value,"T":0.,"A_flux":1-value,"A_vol":1-value,
                       "order_R":value,"order_T":0.,"reference_pseudo_reflection":0.}}


class LedgerTests(unittest.TestCase):
    def test_union_not_sum_and_rank_hours_include_overlap_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            ledger=SolverLedger(Path(directory)/"ledger.json")
            ledger.add_interval(100,200,4,"a","QUALIFIED")
            ledger.add_interval(150,300,4,"b","TIME_LIMIT_UNQUALIFIED")
            ledger.add_interval(350,400,2,"c","CRASHED")
            self.assertEqual(ledger.total_active_seconds(),250)
            self.assertAlmostEqual(ledger.rank_hours(),1100/3600)
            self.assertEqual(ledger.as_dict()["calendar_span_seconds"],300)
            self.assertEqual(SolverLedger(ledger.path).total_active_seconds(),250)

    def test_open_interval_accumulates_until_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            ledger=SolverLedger(Path(directory)/"ledger.json")
            record=ledger.add_interval(100,None,4,"a","RUNNING")
            self.assertEqual(ledger.total_active_seconds(now=250),150)
            ledger.close_interval(record,220,"MEMORY_LIMIT")
            self.assertEqual(ledger.total_active_seconds(now=250),120)

    def test_corrupt_ledger_does_not_reset_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"ledger.json"
            path.write_text('{"version":99,"intervals":[]}')
            with self.assertRaises(ValueError): SolverLedger(path)

    def test_concurrent_ledger_writers_preserve_every_attempt(self):
        with tempfile.TemporaryDirectory() as directory:
            ledger=SolverLedger(Path(directory)/"ledger.json")
            def record(i):
                index=ledger.add_interval(100+i,None,4,str(i),"RUNNING")
                ledger.close_interval(index,200+i,"TIME_LIMIT_UNQUALIFIED")
            with ThreadPoolExecutor(max_workers=2) as pool:
                list(pool.map(record,range(30)))
            loaded=SolverLedger(ledger.path)
            self.assertEqual(len(loaded.data["intervals"]),30)
            self.assertEqual(loaded.total_active_seconds(),129)
            self.assertAlmostEqual(loaded.rank_hours(),12000/3600)


class EvidenceTests(unittest.TestCase):
    def test_complete_qualified_only_and_corruption_detection(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            atomic_json(root/"result.json",qualified())
            (root/"field.bin").write_bytes(b"evidence")
            self.assertTrue(seal_result(root,"implementation-and-case"))
            self.assertIsNotNone(reusable_result(root,"implementation-and-case"))
            self.assertIsNone(reusable_result(root,"different-config"))
            (root/"field.bin").write_bytes(b"tampered")
            self.assertIsNone(reusable_result(root,"implementation-and-case"))

    def test_timeout_is_not_qualified(self):
        with tempfile.TemporaryDirectory() as directory:
            atomic_json(Path(directory)/"result.json",{"status":"TIME_LIMIT_UNQUALIFIED"})
            self.assertFalse(seal_result(directory,"same"))
            self.assertIsNone(reusable_result(directory,"same"))

    def test_cannot_label_missing_stop_qualified(self):
        data=qualified();data["stop"]["converged"]=False
        with self.assertRaises(ValueError): validate_result(data)

    def test_nonfinite_qualified_metric_rejected(self):
        data=qualified();data["metrics"]["A_vol"]=float("nan")
        with self.assertRaises(ValueError): validate_result(data)

    def test_timeout_or_bad_physics_cannot_be_sealed_by_status(self):
        data=qualified();data['stop']['stop_reason']='TIME_LIMIT_UNQUALIFIED'
        with self.assertRaises(ValueError): validate_result(data)
        data=qualified();data['metrics']['A_vol']=.1
        with self.assertRaises(ValueError): validate_result(data)

    def test_case_fingerprint_ignores_identity_not_physics(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            material=root/"material.json";atomic_json(material,{"epsilon":1})
            a=root/"a.json";b=root/"b.json"
            atomic_json(a,{"id":"validation","material_path":str(material),"resolution":20})
            atomic_json(b,{"id":"pilot","material_path":str(material),"resolution":20})
            self.assertEqual(case_fingerprint(a,"code"),case_fingerprint(b,"code"))
            atomic_json(b,{"id":"pilot","material_path":str(material),"resolution":24})
            self.assertNotEqual(case_fingerprint(a,"code"),case_fingerprint(b,"code"))

    def test_lock_excludes_second_writer(self):
        with tempfile.TemporaryDirectory() as directory:
            with FileLock(Path(directory)/"task.lock"):
                with self.assertRaises(RuntimeError):
                    with FileLock(Path(directory)/"task.lock"): pass

    def test_creator_manifest_and_case_hashes_are_enforced(self):
        from ti2d.common import digest_object
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            case={"id":"a","implementation_sha256":"code","resolution":20}
            atomic_json(root/"case.json",case)
            manifest={"implementation_sha256":"code","tasks":[{"id":"a","case_path":"case.json",
                      "case_fingerprint":digest_object({k:v for k,v in case.items() if k!="id"})}]}
            manifest["manifest_sha256"]=digest_object(manifest)
            verify_manifest(manifest,root,"code")
            case["resolution"]=24;atomic_json(root/"case.json",case)
            with self.assertRaisesRegex(ValueError,"CASE_PAYLOAD"):
                verify_manifest(manifest,root,"code")
            manifest["tasks"][0]["id"]="changed"
            with self.assertRaisesRegex(ValueError,"MANIFEST_PAYLOAD"):
                verify_manifest(manifest,root,"code")

    def test_report_does_not_trust_unsealed_qualified_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            taskdir=root/"tasks/a/attempt_1";taskdir.mkdir(parents=True)
            atomic_json(root/"case.json",{"geometry":{}})
            atomic_json(root/"manifest.json",{"tasks":[{"id":"a","phase":"validation","case_path":"case.json"}]})
            atomic_json(taskdir/"result.json",qualified())
            _,_,rows,results=collect(root/"manifest.json")
            self.assertEqual(rows[0]["status"],"UNVERIFIED_EVIDENCE")
            self.assertEqual(results["a"]["status"],"UNVERIFIED_EVIDENCE")

    def test_new_runs_default_to_shared_stage_budget(self):
        manifest={"stage_root":"/stage"}
        self.assertEqual(ledger_path_for_manifest(manifest,"/stage/runs/one"),Path("/stage/budget_ledger.json"))
        self.assertEqual(ledger_path_for_manifest(manifest,"/stage/runs/two"),Path("/stage/budget_ledger.json"))


class ResourceTests(unittest.TestCase):
    def test_container_limits_not_host_memory(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/"cpu.max").write_text("2000000 100000")
            (root/"memory.max").write_text("90000000000")
            (root/"memory.current").write_text("10000000000")
            info=cgroup_resources(root)
            self.assertEqual(info["cpu_quota"],20)
            self.assertEqual(info["memory_available_bytes"],80000000000)

    def test_unbounded_container_memory_refuses_admission(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/"cpu.max").write_text("2000000 100000")
            (root/"memory.max").write_text("max")
            with self.assertRaises(RuntimeError): cgroup_resources(root)

    def test_full_child_tree_memory(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for pid,parent,pages in ((1,0,10),(2,1,20),(3,2,30),(4,0,400)):
                p=root/str(pid);p.mkdir()
                (p/"stat").write_text(f"{pid} (name with spaces) S {parent} 0 0")
                (p/"statm").write_text(f"1000 {pages}")
            with patch("os.sysconf",return_value=4096):
                self.assertEqual(process_tree_rss(1,root),60*4096)

    def test_pair_memory_sum_never_exceeds_group_seventy_percent(self):
        result=pair_admission({"cpu_quota":20,"memory_available_bytes":10*1024**3},4)
        self.assertTrue(result["admitted"])
        self.assertLessEqual(2*result["memory_budget_bytes"],int(.7*10*1024**3))
        self.assertFalse(pair_admission({"cpu_quota":4,"memory_available_bytes":10*1024**3},4)["admitted"])


class ConcurrencyTests(unittest.TestCase):
    def test_slowdown_and_cache_mismatch_return_to_serial(self):
        serial={**qualified(),"reference_reused":False,"queue_process":{"elapsed_s":100}}
        concurrent={**qualified(),"reference_reused":False,"queue_process":{"elapsed_s":131}}
        report=concurrency_assessment(serial,concurrent)
        self.assertTrue(report["compared"])
        self.assertIn("exceeds 30%",report["reason"])
        concurrent["reference_reused"]=True
        self.assertFalse(concurrency_assessment(serial,concurrent)["compared"])

    def test_pair_uses_exactly_two_needed_tasks_with_overlapping_ledger(self):
        with tempfile.TemporaryDirectory() as directory:
            queue=Queue.__new__(Queue)
            queue.run_dir=Path(directory)
            queue.ledger=SolverLedger(Path(directory)/"ledger.json")
            queue.cancelled=False
            barrier=threading.Barrier(2)
            def fake_execute(task,attempt,admission):
                index=queue.ledger.add_interval(time.time(),None,4,task["id"],"RUNNING",attempt)
                barrier.wait(timeout=2)
                time.sleep(.02)
                queue.ledger.close_interval(index,time.time(),"QUALIFIED")
                return {**qualified(),"reference_reused":False,"queue_process":{"elapsed_s":110}},Path(directory)/task["id"]
            queue.execute=fake_execute
            paired=[{"id":"needed_p_minus"},{"id":"needed_s_plus"}]
            serial={**qualified(),"reference_reused":False,"queue_process":{"elapsed_s":100}}
            outputs=queue.execute_pair(paired,serial,{"admitted":True,"memory_budget_bytes":1024**3})
            self.assertEqual(set(outputs),{"needed_p_minus","needed_s_plus"})
            self.assertEqual(len(queue.ledger.data["intervals"]),2)
            self.assertLess(queue.ledger.total_active_seconds(),sum(b-a for a,b in queue.ledger.intervals()))
            self.assertFalse(queue.cancelled)

    def test_failed_probe_stops_peer_without_launching_more(self):
        with tempfile.TemporaryDirectory() as directory:
            queue=Queue.__new__(Queue);queue.run_dir=Path(directory);queue.cancelled=False
            barrier=threading.Barrier(2)
            def fake_execute(task,attempt,admission):
                barrier.wait(timeout=2)
                if task["id"]=="bad":
                    return {"status":"NUMERICAL_FAILED"},Path(directory)/"bad"
                deadline=time.time()+2
                while not queue.cancelled and time.time()<deadline: time.sleep(.001)
                return {"status":"TIME_LIMIT_UNQUALIFIED"},Path(directory)/"peer"
            queue.execute=fake_execute
            outputs=queue.execute_pair([{"id":"bad"},{"id":"peer"}],qualified(),{"admitted":True})
            self.assertTrue(queue.cancelled)
            self.assertEqual(len(outputs),2)


class GateTests(unittest.TestCase):
    def manifest(self):
        return {"production_gate_task_ids":["a","b"],"gate_comparisons":[
            {"kind":kind,"ids":["a","b"],"tolerance":.01,"metrics":["R","T","A_flux"]}
            for kind in ("mesh","mirror","strict")]}

    def test_all_gates_pass_only_with_qualified_evidence(self):
        self.assertTrue(aggregate_gate(self.manifest(),{"a":qualified(.2),"b":qualified(.201)})["passed"])
        self.assertFalse(aggregate_gate(self.manifest(),{"a":qualified(.2)})["passed"])

    def test_exact_threshold_fails(self):
        manifest=self.manifest();manifest["gate_comparisons"][0]["tolerance"]=.125
        self.assertFalse(aggregate_gate(manifest,{"a":qualified(.25),"b":qualified(.375)})["passed"])

    def test_missing_comparison_cannot_unlock_pilot(self):
        manifest=self.manifest();manifest["gate_comparisons"].pop()
        self.assertFalse(aggregate_gate(manifest,{"a":qualified(),"b":qualified()})["passed"])

    def test_unknown_bound_does_not_invent_contrast(self):
        rows=[{"phase":"pilot","geometry_id":"baseline","polarization":"p","emission_theta_deg":theta,
               "A_flux":a,"A_error_bound":None,"status":"QUALIFIED"} for theta,a in ((-30,.2),(30,.3))]
        pair=[p for p in paired_metrics(rows,True) if p["polarization"]=="p"][0]
        self.assertIsNone(pair["normalized_contrast"])
        self.assertAlmostEqual(pair["angle_difference"],.1)

    def test_contrast_floor_and_unresolved_effect(self):
        rows=[{"phase":"pilot","geometry_id":"baseline","polarization":"p","emission_theta_deg":theta,
               "A_flux":a,"A_error_bound":.01,"status":"QUALIFIED"} for theta,a in ((-30,.02),(30,.03))]
        pair=[p for p in paired_metrics(rows,True) if p["polarization"]=="p"][0]
        self.assertIsNone(pair["normalized_contrast"])
        self.assertIn("cannot resolve",pair["interpretation"])


if __name__ == "__main__":
    unittest.main()
