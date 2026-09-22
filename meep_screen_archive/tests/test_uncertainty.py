"""Pure saved-evidence tests; no Meep or external files are needed."""

from copy import deepcopy
import json
import unittest

from ti2d.uncertainty import derive_error_evidence


def fixture(case="structure"):
    baseline = {
        "id": "base", "geometry_id": "baseline", "geometry": {
            "period_um": 50.0, "surface_fill": 0.692820323, "axis_length_um": 40.0,
            "tilt_deg": 30.0, "thickness_um": 60.0},
        "case": case, "material_sha256": "a" * 64, "implementation_sha256": "b" * 64,
        "wavelength_um": 10.5, "polarization": "p", "emission_theta_deg": 30.0,
        "resolution": 20, "courant": 0.25, "strict_stop": False, "intensity_decay": 1e-4,
    }
    tasks, results = [], {}
    for name, resolution, strict, metrics in (
        ("mesh16", 16, False, {"R": 0.31, "T": 0.1, "A_flux": 0.59, "A_vol": 0.591}),
        ("base", 20, False, {"R": 0.30, "T": 0.1, "A_flux": 0.60, "A_vol": 0.602}),
        ("mesh24", 24, False, {"R": 0.295, "T": 0.1, "A_flux": 0.605, "A_vol": 0.604}),
        ("strict", 20, True, {"R": 0.298, "T": 0.099, "A_flux": 0.603, "A_vol": 0.603}),
    ):
        config = deepcopy(baseline)
        config.update(id=name, resolution=resolution, strict_stop=strict,
                      intensity_decay=1e-6 if strict else 1e-4)
        task = {"id": name, "phase": "validation", "geometry_id": "baseline", "resolution": resolution,
                "polarization": "p", "emission_theta_deg": 30.0, "strict_stop": strict}
        tasks.append(task)
        results[name] = {
            "status": "QUALIFIED", "case_config": config, "metrics": metrics,
            "stop": {"converged": True,
                     "windows": [{"start": start, "end": start + 20} for start in (0, 20, 40)],
                     "trace": [{"time": time, "R": 0.30 + time * 1e-5, "T": 0.10 + time * 2e-6}
                               for time in (10, 20, 30, 40, 50, 60)]},
        }
        if case == "flat":
            results[name]["analytical"] = {"R": 0.301, "T": 0.101, "A": 0.598}
    manifest = {"tasks": tasks, "material_sha256": "a" * 64, "implementation_sha256": "b" * 64,
                "gate_comparisons": [
                    {"kind": "mesh", "ids": ["mesh16", "base", "mesh24"], "metrics": ["R", "T", "A_flux"]},
                    {"kind": "strict", "ids": ["base", "strict"], "metrics": ["R", "T", "A_flux"]},
                ]}
    return manifest, results


def rows_by_id(manifest, results):
    return {row["id"]: row for row in derive_error_evidence(manifest, results)}


class ErrorEvidenceTests(unittest.TestCase):
    def test_exact_grid_strict_absorption_and_tail_sum_is_empirical_only(self):
        manifest, results = fixture()
        rows = rows_by_id(manifest, results)
        base = rows["base"]
        self.assertAlmostEqual(base["mesh_span_R"], 0.015)
        self.assertEqual(base["mesh_span_T"], 0)
        self.assertAlmostEqual(base["mesh_span_A_flux"], 0.015)
        self.assertAlmostEqual(base["strict_stop_delta_A_flux"], 0.003)
        self.assertAlmostEqual(base["absorption_delta_A"], 0.002)
        self.assertAlmostEqual(base["tail_R_range"], 0.0005)
        self.assertAlmostEqual(base["tail_T_range"], 0.0001)
        self.assertAlmostEqual(base["A_error_estimate"], 0.0206)
        self.assertEqual(base["missing_coverage"], [])
        self.assertEqual(base["estimate_kind"], "empirical_not_rigorous")
        self.assertFalse(base["validated"])
        self.assertFalse(base["formal_error_bound_available"])
        self.assertTrue(base["identity_verified"])
        self.assertEqual(base["A_error_estimate"], sum(base["estimate_components"].values()))
        json.dumps(list(rows.values()), allow_nan=False)

    def test_different_resolution_never_inherits_strict_stop_evidence(self):
        manifest, results = fixture()
        rows = rows_by_id(manifest, results)
        for name in ("mesh16", "mesh24"):
            self.assertAlmostEqual(rows[name]["mesh_span_A_flux"], 0.015)
            self.assertIsNone(rows[name]["strict_stop_delta_A_flux"])
            self.assertIsNone(rows[name]["A_error_estimate"])
            self.assertIn("strict_stop", rows[name]["missing_coverage"])

    def test_strict_variant_never_inherits_mesh_at_different_stopping_controls(self):
        manifest, results = fixture()
        strict = rows_by_id(manifest, results)["strict"]
        self.assertAlmostEqual(strict["strict_stop_delta_A_flux"], 0.003)
        self.assertIsNone(strict["mesh_span_A_flux"])
        self.assertIsNone(strict["A_error_estimate"])
        self.assertIn("mesh", strict["missing_coverage"])

    def test_exact_reused_alias_has_same_evidence_with_own_task_identity(self):
        manifest, results = fixture()
        alias = {**manifest["tasks"][1], "id": "pilot_baseline", "phase": "pilot", "reused_from": "base"}
        manifest["tasks"].append(alias)
        results["pilot_baseline"] = results["base"]
        rows = rows_by_id(manifest, results)
        self.assertEqual(rows["pilot_baseline"]["A_error_estimate"], rows["base"]["A_error_estimate"])
        self.assertEqual(rows["pilot_baseline"]["id"], "pilot_baseline")
        self.assertEqual(rows["pilot_baseline"]["phase"], "pilot")

    def test_mesh_boundary_comparison_is_supported(self):
        manifest, results = fixture()
        manifest["gate_comparisons"][0] = {"kind": "mesh_boundary", "ids": ["base", "mesh24"]}
        base = rows_by_id(manifest, results)["base"]
        self.assertAlmostEqual(base["mesh_span_A_flux"], 0.005)
        self.assertEqual(base["mesh_comparisons"][0]["kind"], "mesh_boundary")

    def test_mirror_comparison_is_not_mesh_or_strict_evidence(self):
        manifest, results = fixture()
        manifest["gate_comparisons"] = [{"kind": "mirror", "ids": ["base", "mesh24"]}]
        base = rows_by_id(manifest, results)["base"]
        self.assertIsNone(base["mesh_span_A_flux"])
        self.assertIsNone(base["strict_stop_delta_A_flux"])
        self.assertIsNone(base["A_error_estimate"])

    def test_no_cross_angle_polarization_wavelength_or_material_transfer(self):
        for field, value in (("emission_theta_deg", -30.0), ("polarization", "s"),
                             ("wavelength_um", 10.6), ("material_sha256", "c" * 64),
                             ("implementation_sha256", "d" * 64)):
            for member, evidence in (("mesh24", "mesh_span_A_flux"), ("strict", "strict_stop_delta_A_flux")):
                manifest, results = fixture()
                results[member]["case_config"][field] = value
                for task in manifest["tasks"]:
                    if task["id"] == member and field in task:
                        task[field] = value
                with self.subTest(field=field, evidence=evidence):
                    base = rows_by_id(manifest, results)["base"]
                    self.assertIsNone(base[evidence])
                    self.assertIsNone(base["A_error_estimate"])

    def test_geometry_content_and_other_numerical_controls_must_match_exactly(self):
        for key, value in (("geometry", {"period_um": 58.0}), ("courant", 0.2),
                           ("case", "flat"), ("source_cutoff", 6.0)):
            manifest, results = fixture()
            results["mesh24"]["case_config"][key] = value
            with self.subTest(key=key):
                base = rows_by_id(manifest, results)["base"]
                self.assertIsNone(base["mesh_span_A_flux"])
                self.assertIsNone(base["A_error_estimate"])

    def test_manifest_label_disagreement_invalidates_identity(self):
        manifest, results = fixture()
        results["base"]["case_config"]["emission_theta_deg"] = -30.0
        base = rows_by_id(manifest, results)["base"]
        self.assertFalse(base["identity_verified"])
        self.assertIn("emission_theta_deg", base["identity_error"])
        self.assertIn("exact_case_identity", base["missing_coverage"])
        self.assertIsNone(base["A_error_estimate"])

    def test_failed_results_are_excluded_and_do_not_contribute_comparison_values(self):
        for status in ("TIME_LIMIT_UNQUALIFIED", "CRASHED", "NUMERICAL_FAILED", "SMOKE_ONLY"):
            manifest, results = fixture()
            results["mesh24"]["status"] = status
            results["mesh24"]["metrics"]["A_flux"] = 999
            with self.subTest(status=status):
                rows = rows_by_id(manifest, results)
                self.assertNotIn("mesh24", rows)
                self.assertIsNone(rows["base"]["mesh_span_A_flux"])
                self.assertIsNone(rows["base"]["A_error_estimate"])

    def test_results_outside_manifest_are_ignored(self):
        manifest, results = fixture()
        results["old_legacy_case"] = deepcopy(results["base"])
        self.assertNotIn("old_legacy_case", rows_by_id(manifest, results))

    def test_missing_or_nonfinite_required_component_prevents_combined_estimate(self):
        for result_id, metric in (("base", "A_vol"), ("mesh24", "A_flux"), ("strict", "A_flux")):
            for bad in (None, float("nan"), float("inf"), True):
                manifest, results = fixture()
                results[result_id]["metrics"][metric] = bad
                with self.subTest(result_id=result_id, metric=metric, bad=bad):
                    rows = rows_by_id(manifest, results)
                    self.assertIsNone(rows["base"]["A_error_estimate"])
                    json.dumps(list(rows.values()), allow_nan=False)

    def test_flat_fresnel_difference_is_added_only_to_exact_flat_result(self):
        manifest, results = fixture(case="flat")
        base = rows_by_id(manifest, results)["base"]
        self.assertAlmostEqual(base["flat_fresnel_delta_R"], 0.001)
        self.assertAlmostEqual(base["flat_fresnel_delta_T"], 0.001)
        self.assertAlmostEqual(base["flat_fresnel_delta_A_flux"], 0.002)
        self.assertAlmostEqual(base["A_error_estimate"], 0.0226)
        manifest, results = fixture()
        results["base"]["analytical"] = {"R": 0, "T": 0, "A": 1}
        base = rows_by_id(manifest, results)["base"]
        self.assertIsNone(base["flat_fresnel_delta_A_flux"])

    def test_missing_optional_fresnel_or_tail_is_explicit_not_an_invented_zero(self):
        manifest, results = fixture(case="flat")
        results["base"].pop("analytical")
        results["base"]["stop"].pop("trace")
        base = rows_by_id(manifest, results)["base"]
        self.assertEqual(set(base["missing_coverage"]), {"flat_fresnel", "tail"})
        self.assertIsNone(base["flat_fresnel_delta_A_flux"])
        self.assertIsNone(base["tail_A_estimate"])
        self.assertAlmostEqual(base["A_error_estimate"], 0.020)
        self.assertNotIn("tail_A_estimate", base["estimate_components"])
        self.assertNotIn("flat_fresnel_delta_A_flux", base["estimate_components"])

    def test_tail_uses_cross_window_range_instead_of_zero_within_window_ranges(self):
        manifest, results = fixture()
        trace = results["base"]["stop"]["trace"]
        for index, sample in enumerate(trace):
            sample["R"] = 0.3 + (index // 2) * 0.001
            sample["T"] = 0.1
        for window in results["base"]["stop"]["windows"]:
            window["R_range"] = window["T_range"] = 0
        base = rows_by_id(manifest, results)["base"]
        self.assertAlmostEqual(base["tail_R_range"], 0.002)
        self.assertAlmostEqual(base["tail_A_estimate"], 0.002)

    def test_tail_uses_exactly_last_three_complete_windows(self):
        manifest, results = fixture()
        stop = results["base"]["stop"]
        stop["windows"].insert(0, {"start": -20, "end": 0})
        stop["trace"].insert(0, {"time": -10, "R": 1000, "T": 1000})
        stop["trace"].append({"time": 62, "R": 1000, "T": 1000})
        base = rows_by_id(manifest, results)["base"]
        self.assertEqual(base["tail_start"], 0)
        self.assertEqual(base["tail_end"], 60)
        self.assertEqual(base["tail_samples"], 6)
        self.assertAlmostEqual(base["tail_R_range"], 0.0005)

    def test_sparse_nonfinite_or_discontinuous_tail_has_missing_coverage(self):
        for modification in ("sparse", "nonfinite", "discontinuous"):
            manifest, results = fixture()
            stop = results["base"]["stop"]
            if modification == "sparse":
                stop["trace"] = stop["trace"][::2]
            elif modification == "nonfinite":
                stop["trace"][2]["R"] = float("nan")
            else:
                stop["windows"][1]["start"] += 1
            with self.subTest(modification=modification):
                base = rows_by_id(manifest, results)["base"]
                self.assertIsNone(base["tail_R_range"])
                self.assertIn("tail", base["missing_coverage"])

    def test_missing_config_cannot_borrow_a_named_geometrys_evidence(self):
        manifest, results = fixture()
        results["base"].pop("case_config")
        base = rows_by_id(manifest, results)["base"]
        self.assertFalse(base["identity_verified"])
        self.assertIsNone(base["A_error_estimate"])
        self.assertIn("exact_case_identity", base["missing_coverage"])
        self.assertAlmostEqual(base["absorption_delta_A"], 0.002)

    def test_declared_strict_pair_must_have_genuinely_tighter_decay(self):
        manifest, results = fixture()
        results["strict"]["case_config"]["intensity_decay"] = 1e-3
        base = rows_by_id(manifest, results)["base"]
        self.assertIsNone(base["strict_stop_delta_A_flux"])
        self.assertIsNone(base["A_error_estimate"])

    def test_unsampled_resolution_does_not_inherit_a_surrounding_mesh_span(self):
        manifest, results = fixture()
        task = {**manifest["tasks"][1], "id": "unsampled22", "resolution": 22}
        manifest["tasks"].append(task)
        results["unsampled22"] = deepcopy(results["base"])
        results["unsampled22"]["case_config"].update(id="unsampled22", resolution=22)
        row = rows_by_id(manifest, results)["unsampled22"]
        self.assertIsNone(row["mesh_span_A_flux"])
        self.assertIsNone(row["A_error_estimate"])

    def test_input_objects_are_not_modified(self):
        manifest, results = fixture()
        before_manifest, before_results = deepcopy(manifest), deepcopy(results)
        derive_error_evidence(manifest, results)
        self.assertEqual(manifest, before_manifest)
        self.assertEqual(results, before_results)


if __name__ == "__main__":
    unittest.main()
