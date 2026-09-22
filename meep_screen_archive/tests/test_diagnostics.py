"""Analytical tests: no Meep import, FDTD execution, or calibration fixtures."""

import json
import math
import unittest

import numpy as np

from ti2d.diagnostics import (
    StopTracker,
    direction_metrics,
    fresnel_slab,
    integrate_absorption,
    periodic_projection,
)


def plane_fields(x, period, kx, f, pol, waves):
    """Synthesize exact Maxwell plane waves independently of the projection."""
    result = {name: np.zeros(x.shape, dtype=complex)
              for name in (("Ex", "Hz") if pol == "p" else ("Ez", "Hx"))}
    for m, upward, downward in waves:
        transverse = kx + m / period
        normal = np.sqrt(complex(f * f - transverse * transverse)) / f
        phase = np.exp(2j * np.pi * transverse * x)
        if pol == "p":
            result["Hz"] += (upward + downward) * phase
            result["Ex"] += normal * (downward - upward) * phase
        else:
            result["Ez"] += (upward + downward) * phase
            result["Hx"] += normal * (upward - downward) * phase
    return result


class FresnelTests(unittest.TestCase):
    def test_vacuum_is_transparent_at_every_angle(self):
        for pol in ("s", "p"):
            for angle in (0, 23, -57, 82):
                with self.subTest(pol=pol, angle=angle):
                    result = fresnel_slab(1, 10, angle, pol)
                    self.assertAlmostEqual(result["R"], 0, places=12)
                    self.assertAlmostEqual(result["T"], 1, places=12)
                    self.assertEqual(result["A"], 0)

    def test_zero_thickness_is_transparent_even_for_metal(self):
        for pol in ("s", "p"):
            result = fresnel_slab(-150 + 39j, 10, 63, pol, thickness_um=0)
            self.assertEqual((result["R"], result["T"], result["A"]), (0, 1, 0))

    def test_quarter_wave_film_matches_independent_formula(self):
        n = 2
        # Air/n/air quarter-wave film: R=((n²-1)/(n²+1))².
        expected = ((n * n - 1) / (n * n + 1)) ** 2
        for pol in ("s", "p"):
            result = fresnel_slab(n * n, 8, 0, pol, thickness_um=1)
            self.assertAlmostEqual(result["R"], expected, places=13)
            self.assertAlmostEqual(result["T"], 1 - expected, places=13)
            self.assertEqual(result["A"], 0)

    def test_lossless_oblique_energy_and_angle_reversal(self):
        for pol in ("s", "p"):
            plus = fresnel_slab(2.7, 8.3, 46, pol, thickness_um=1.9)
            minus = fresnel_slab(2.7, 8.3, -46, pol, thickness_um=1.9)
            self.assertAlmostEqual(plus["R"] + plus["T"], 1, places=12)
            self.assertEqual(plus["R"], minus["R"])

    def test_thick_metal_is_finite_and_matches_single_interface(self):
        eps = -200 + 70j
        angle = math.radians(35)
        q0 = math.cos(angle)
        q1 = np.sqrt(eps - math.sin(angle) ** 2)
        for pol in ("s", "p"):
            result = fresnel_slab(eps, 0.1, 35, pol, thickness_um=1e6)
            r = (q0 - q1) / (q0 + q1) if pol == "s" else (eps * q0 - q1) / (eps * q0 + q1)
            self.assertAlmostEqual(result["R"], abs(r) ** 2, places=12)
            self.assertEqual(result["T"], 0)
            self.assertGreater(result["A"], 0)
            self.assertTrue(np.isfinite(result["A"]))

    def test_finite_metal_transmission_and_loss(self):
        for pol in ("s", "p"):
            result = fresnel_slab(-4 + 2j, 10, 32, pol, thickness_um=0.02)
            self.assertGreater(result["T"], 0)
            self.assertGreater(result["A"], 0)
            self.assertAlmostEqual(result["R"] + result["T"] + result["A"], 1)

    def test_finite_film_matches_independent_boundary_equation_solution(self):
        # Solve four continuity equations for r, internal A/B and transmitted t.
        # This independent construction does not use a Fresnel recursion.
        eps, wavelength, theta, thickness = -3 + 2j, 9.3, 41, 0.47
        q0 = math.cos(math.radians(theta))
        q1 = np.sqrt(eps - math.sin(math.radians(theta)) ** 2)
        phase = np.exp(2j * math.pi * q1 * thickness / wavelength)
        for pol in ("s", "p"):
            normal_admittance = q1 if pol == "s" else q1 / eps
            matrix = np.array([
                [1, -1, -1, 0],
                [-q0, -normal_admittance, normal_admittance, 0],
                [0, phase, 1 / phase, -1],
                [0, normal_admittance * phase, -normal_admittance / phase, -q0],
            ], dtype=complex)
            solution = np.linalg.solve(matrix, [-1, -q0, 0, 0])
            result = fresnel_slab(eps, wavelength, theta, pol, thickness_um=thickness)
            self.assertAlmostEqual(result["R"], abs(solution[0]) ** 2, places=12)
            self.assertAlmostEqual(result["T"], abs(solution[3]) ** 2, places=12)

    def test_zero_normal_wavevector_limits_are_finite(self):
        for angle in (0, 30):
            eps = math.sin(math.radians(angle)) ** 2
            for pol in ("s", "p"):
                exact = fresnel_slab(eps, 8, angle, pol, thickness_um=1)
                nearby = fresnel_slab(eps + 1e-8, 8, angle, pol, thickness_um=1)
                self.assertAlmostEqual(exact["R"] + exact["T"], 1, places=10)
                self.assertAlmostEqual(exact["R"], nearby["R"], places=6)

    def test_invalid_inputs(self):
        for kwargs in ({"eps": 1 - 1j}, {"lambda_um": 0}, {"theta_deg": 90},
                       {"pol": "TE"}, {"thickness_um": -1}, {"eps": np.nan},
                       {"theta_deg": np.complex128(1 + 1j)}):
            arguments = dict(eps=2 + 1j, lambda_um=10, theta_deg=0, pol="s")
            arguments.update(kwargs)
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                fresnel_slab(**arguments)


class ProjectionTests(unittest.TestCase):
    def setUp(self):
        self.period = 20.0
        self.f = 0.21
        self.x = -10 + (np.arange(64) + 0.5) * self.period / 64

    def project(self, pol, waves, kx=0, x=None):
        x = self.x if x is None else x
        fields = plane_fields(x, self.period, kx, self.f, pol, waves)
        return periodic_projection(x, fields, self.period, kx, self.f, pol)

    def test_normal_incidence_both_polarizations_and_directions(self):
        for pol in ("s", "p"):
            for upward, downward in ((1 + 0.3j, 0), (0, 2 - 0.2j)):
                with self.subTest(pol=pol, upward=upward):
                    result = self.project(pol, [(0, upward, downward)])
                    self.assertAlmostEqual(result["up_power"], self.period * abs(upward) ** 2, places=11)
                    self.assertAlmostEqual(result["down_power"], self.period * abs(downward) ** 2, places=11)
                    self.assertAlmostEqual(result["direct_net_up_power"], result["net_up_power"], places=11)
                    zero = next(order for order in result["orders"] if order["m"] == 0)
                    self.assertAlmostEqual(abs(zero["up_amplitude"] - upward), 0, places=12)
                    self.assertAlmostEqual(abs(zero["down_amplitude"] - downward), 0, places=12)

    def test_oblique_bloch_phase_and_negative_incident_kx(self):
        kx = -self.f * math.sin(math.radians(42))
        q = math.cos(math.radians(42))
        for pol in ("s", "p"):
            result = self.project(pol, [(0, 0, 1.7j)], kx=kx)
            self.assertAlmostEqual(result["down_power"], self.period * q * 1.7 ** 2, places=11)
            self.assertLess(result["up_power"], 1e-27)
            zero = next(order for order in result["orders"] if order["m"] == 0)
            self.assertAlmostEqual(zero["kx"], kx)

    def test_multiple_orders_and_both_directions(self):
        kx = -0.031
        waves = [(-2, 0.4 + 0.9j, 0.1), (0, 0.2j, 1.0), (3, 0.8, -0.1j)]
        expected_up = sum(self.period * math.sqrt(1 - ((kx + m / self.period) / self.f) ** 2) * abs(up) ** 2
                          for m, up, down in waves)
        expected_down = sum(self.period * math.sqrt(1 - ((kx + m / self.period) / self.f) ** 2) * abs(down) ** 2
                            for m, up, down in waves)
        for pol in ("s", "p"):
            result = self.project(pol, waves, kx=kx)
            self.assertAlmostEqual(result["up_power"], expected_up, places=11)
            self.assertAlmostEqual(result["down_power"], expected_down, places=11)
            self.assertAlmostEqual(result["net_up_power"], result["direct_net_up_power"], places=11)

    def test_reference_subtraction_must_be_complex_before_projection(self):
        kx = -0.061
        for pol in ("s", "p"):
            incident = plane_fields(self.x, self.period, kx, self.f, pol, [(0, 0, 2 + 0.4j)])
            total = plane_fields(self.x, self.period, kx, self.f, pol, [(0, 0.1 - 0.6j, 2 + 0.4j)])
            scattered = {name: total[name] - incident[name] for name in incident}
            result = periodic_projection(self.x, scattered, self.period, kx, self.f, pol)
            expected = self.period * math.sqrt(1 - (kx / self.f) ** 2) * abs(0.1 - 0.6j) ** 2
            self.assertAlmostEqual(result["up_power"], expected, places=11)
            self.assertLess(result["down_power"], 1e-28)

    def test_signed_direct_power_retains_downward_sign(self):
        for pol in ("s", "p"):
            result = self.project(pol, [(0, 0, 1)])
            self.assertLess(result["net_up_power"], 0)
            self.assertLess(result["direct_net_up_power"], 0)
            self.assertGreater(result["down_power"], 0)

    def test_endpoint_grid_not_double_counted_and_bloch_preserved(self):
        x = np.linspace(-10, 10, 65)
        for pol in ("s", "p"):
            endpoint = self.project(pol, [(0, 1j, 0), (2, 0.3, 0)], kx=-0.021, x=x)
            midpoint = self.project(pol, [(0, 1j, 0), (2, 0.3, 0)], kx=-0.021)
            self.assertTrue(endpoint["duplicate_endpoint_removed"])
            self.assertEqual(endpoint["samples"], 64)
            self.assertAlmostEqual(endpoint["up_power"], midpoint["up_power"], places=11)

    def test_origin_translation_preserves_power(self):
        for pol in ("s", "p"):
            a = self.project(pol, [(1, 1, 0.4j)], kx=-0.04)
            b = self.project(pol, [(1, 1, 0.4j)], kx=-0.04, x=self.x + 113.7)
            self.assertAlmostEqual(a["up_power"], b["up_power"], places=11)

    def test_evanescent_order_has_zero_outgoing_power(self):
        for pol in ("s", "p"):
            result = self.project(pol, [(7, 1, 0)])
            evanescent = next(order for order in result["orders"] if order["m"] == 7)
            self.assertFalse(evanescent["propagating"])
            self.assertGreater(evanescent["q_imag"], 0)
            self.assertEqual(evanescent["up_power"], 0)
            self.assertLess(result["up_power"], 1e-27)

    def test_exact_cutoff_is_flagged_without_division(self):
        fields = {"Ez": np.zeros(64), "Hx": np.zeros(64)}
        result = periodic_projection(self.x, fields, self.period, 0, 0.2, "s")
        self.assertEqual(result["cutoff_orders"], [-4, 4])
        for order in result["orders"]:
            if order["cutoff"]:
                self.assertIsNone(order["up_amplitude"])
                self.assertEqual(order["up_power"], 0)

    def test_invalid_grid_and_field_arrays(self):
        good = plane_fields(self.x, self.period, 0, self.f, "s", [(0, 1, 0)])
        cases = [
            (self.x[:-1], good),
            (self.x, {"Ez": good["Ez"][:, None], "Hx": good["Hx"][:, None]}),
            (self.x, {"Ez": good["Ez"], "Hx": good["Hx"], "extra": np.zeros((64, 1))}),
            (self.x, {"Ez": good["Ez"], "Hx": np.full(64, np.nan)}),
            (self.x * 0.7, good),
            (np.where(np.arange(64) == 3, self.x + 0.1, self.x), good),
        ]
        for x, fields in cases:
            with self.subTest(shape=x.shape), self.assertRaises(ValueError):
                periodic_projection(x, fields, self.period, 0, self.f, "s")

    def test_endpoint_must_match_bloch_phase(self):
        x = np.linspace(-10, 10, 65)
        fields = plane_fields(x, self.period, 0.03, self.f, "s", [(0, 1, 0)])
        fields["Ez"][-1] += 0.01
        with self.assertRaisesRegex(ValueError, "Bloch"):
            periodic_projection(x, fields, self.period, 0.03, self.f, "s")

    def test_underresolved_grid_is_rejected(self):
        x = np.array([-5.0, 5.0])
        fields = {"Ez": np.ones(2), "Hx": np.ones(2)}
        with self.assertRaisesRegex(ValueError, "resolve"):
            periodic_projection(x, fields, self.period, 0, self.f, "s")


class AbsorptionTests(unittest.TestCase):
    def pair(self):
        E = np.array([[1 + 1j, 2], [0.5j, 3 - 1j]])
        return {"name": "Ex", "E": E, "D": (2 + 0.4j) * E,
                "weights": np.full((2, 2), 0.3),
                "ti_mask": np.array([[True, False], [True, True]]),
                "coords": {"x": np.array([-0.5, 0.5]), "y": np.array([0, 0.3])}}

    def test_known_uniform_dielectric_loss_and_material_split(self):
        pair = self.pair()
        f, incident = 0.12, 3.7
        result = integrate_absorption([pair], f, incident)
        density = 2 * math.pi * f * 0.4 * np.abs(pair["E"]) ** 2
        expected = float(np.sum(density * pair["weights"]))
        ti = float(np.sum((density * pair["weights"])[pair["ti_mask"]]))
        self.assertAlmostEqual(result["absorbed_power"], expected, places=13)
        self.assertAlmostEqual(result["A"], expected / incident, places=13)
        self.assertAlmostEqual(result["ti_absorbed_power"], ti, places=13)
        self.assertAlmostEqual(result["A_ti"] + result["A_outside_ti"], result["A"], places=13)
        self.assertFalse(result["empirical_calibration"])
        json.dumps(result, allow_nan=False)

    def test_multiple_components_add_without_calibration(self):
        first = self.pair()
        second = self.pair()
        second["name"] = "Ey"
        one = integrate_absorption([first], 0.1, 2)
        both = integrate_absorption([first, second], 0.1, 2)
        self.assertAlmostEqual(both["A"], 2 * one["A"])

    def test_negative_loss_is_preserved_as_diagnostic(self):
        pair = self.pair()
        pair["D"] = (2 - 0.4j) * pair["E"]
        result = integrate_absorption([pair], 0.1, 1)
        self.assertLess(result["A"], 0)
        self.assertLess(result["components"][0]["negative_absorbed_power"], 0)
        self.assertEqual(result["components"][0]["positive_absorbed_power"], 0)

    def test_lossless_vacuum_absorption_zero(self):
        pair = self.pair()
        pair["D"] = pair["E"].copy()
        self.assertEqual(integrate_absorption([pair], 0.1, 1)["A"], 0)

    def test_complex_reference_phase_does_not_change_absorption(self):
        before = integrate_absorption([self.pair()], 0.1, 1)
        after_pair = self.pair()
        for name in ("E", "D"):
            after_pair[name] *= np.exp(0.43j)
        after = integrate_absorption([after_pair], 0.1, 1)
        self.assertAlmostEqual(before["A"], after["A"], places=13)

    def test_no_shape_broadcasting_or_squeezing(self):
        for key in ("D", "weights", "ti_mask"):
            for shape in ((2, 1), (2, 2, 1), (4,)):
                pair = self.pair()
                pair[key] = np.ones(shape)
                with self.subTest(key=key, shape=shape), self.assertRaisesRegex(ValueError, "shapes"):
                    integrate_absorption([pair], 0.1, 1)

    def test_invalid_weights_masks_coordinates_and_values(self):
        mutations = [
            ("weights", np.full((2, 2), -1.0)),
            ("weights", np.zeros((2, 2))),
            ("weights", np.ones((2, 2), dtype=complex)),
            ("ti_mask", np.full((2, 2), 0.5)),
            ("E", np.full((2, 2), np.nan)),
            ("D", np.full((2, 2), np.inf)),
            ("coords", {"x": np.arange(3), "y": np.arange(2)}),
            ("coords", {"x": np.arange(2)}),
        ]
        for key, value in mutations:
            pair = self.pair()
            pair[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                integrate_absorption([pair], 0.1, 1)
        for incident in (0, -1, np.nan):
            with self.assertRaises(ValueError):
                integrate_absorption([self.pair()], 0.1, incident)


class StopTrackerTests(unittest.TestCase):
    def start(self, **kwargs):
        tracker = StopTracker(source_end=10, **kwargs)
        tracker.sample(0, {"slot": 4, "air": 1}, 9, R=0.3, T=0.1)
        tracker.sample(10, {"slot": 1, "air": 0.1}, 1, R=0.3, T=0.1)
        return tracker

    def sample_tail(self, tracker, start=12, end=70, power=1e-5, envelope=1e-5, R=0.3, T=0.1):
        for time in range(start, end + 1, 2):
            tracker.sample(time, {"slot": power, "air": power}, envelope, R=R, T=T)
        return tracker

    def test_three_complete_windows_required(self):
        tracker = self.start()
        self.sample_tail(tracker, end=68)
        self.assertFalse(tracker.converged)
        self.assertTrue(tracker.sample(70, {"slot": 1e-5, "air": 1e-5}, 1e-5, R=0.3, T=0.1))
        self.assertEqual(tracker.summary()["consecutive_passed_windows"], 3)
        self.assertFalse(tracker.summary()["timeout_is_convergence"])
        json.dumps(tracker.summary(), allow_nan=False)
        json.dumps(tracker.trace, allow_nan=False)

    def test_timeout_does_not_make_persistent_fields_converge(self):
        tracker = self.start()
        self.sample_tail(tracker, end=250, power=0.2, envelope=0.4)
        self.assertFalse(tracker.converged)
        self.assertEqual(tracker.summary()["status"], "not_converged")

    def test_non_pml_envelope_catches_mode_missed_by_probes(self):
        tracker = self.start()
        self.sample_tail(tracker, end=170, power=1e-10, envelope=0.2)
        self.assertFalse(tracker.converged)
        self.assertLess(tracker.trace[-1]["probe_ratio"]["slot"], tracker.decay_threshold)
        self.assertGreater(tracker.trace[-1]["envelope_ratio"], tracker.decay_threshold)

    def test_one_undecayed_probe_blocks_stopping(self):
        tracker = self.start()
        for time in range(12, 92, 2):
            tracker.sample(time, {"slot": 1e-8, "air": 0.1}, 1e-8, R=0.3, T=0.1)
        self.assertFalse(tracker.converged)

    def test_historical_peak_not_tail_relative(self):
        tracker = self.start()
        self.sample_tail(tracker, power=5e-5, envelope=5e-5)
        self.assertTrue(tracker.converged)
        self.assertEqual(tracker.trace[-1]["probe_peak"]["slot"], 4)
        # The excitation envelope peak supplies the default reference floor.
        self.assertEqual(tracker.trace[-1]["probe_normalization_denominator"]["slot"], 9)
        self.assertAlmostEqual(tracker.trace[-1]["probe_ratio"]["slot"], 5e-5 / 9)

    def test_probe_above_reference_floor_keeps_own_historical_peak(self):
        tracker = self.start(reference_intensity_peak=1,
                             reference_provenance="qualified empty-cell source-on envelope")
        self.sample_tail(tracker, power=5e-5, envelope=5e-5)
        self.assertTrue(tracker.converged)
        self.assertEqual(tracker.trace[-1]["probe_normalization_denominator"]["slot"], 4)
        self.assertAlmostEqual(tracker.trace[-1]["probe_ratio"]["slot"], 5e-5 / 4)

    def test_absolute_rt_stability_at_tiny_transmission(self):
        tracker = self.start()
        for time in range(12, 72, 2):
            tracker.sample(time, {"slot": 1e-5, "air": 1e-5}, 1e-5,
                           R=0.3, T=(1e-9 if time % 4 else 1e-5))
        self.assertTrue(tracker.converged)

    def test_rt_oscillation_fails_even_when_fields_decay(self):
        tracker = self.start()
        for time in range(12, 112, 2):
            tracker.sample(time, {"slot": 1e-8, "air": 1e-8}, 1e-8,
                           R=0.3 + (0.003 if time % 4 else 0), T=0.1)
        self.assertFalse(tracker.converged)

    def test_accumulated_rt_drift_across_three_windows_fails(self):
        tracker = self.start()
        for time in range(12, 72, 2):
            tracker.sample(time, {"slot": 1e-8, "air": 1e-8}, 1e-8,
                           R=0.3 + time * 3e-5, T=0.1)
        self.assertTrue(all(window["passed"] for window in tracker.summary()["windows"]))
        self.assertFalse(tracker.converged)

    def test_missing_rt_never_passes(self):
        tracker = self.start()
        self.sample_tail(tracker, end=150, R=None, T=None)
        self.assertFalse(tracker.converged)

    def test_sparse_samples_do_not_create_three_valid_windows(self):
        tracker = self.start()
        for time in (30, 50, 70, 90):
            tracker.sample(time, {"slot": 1e-8, "air": 1e-8}, 1e-8, R=0.3, T=0.1)
        self.assertFalse(tracker.converged)

    def test_identically_zero_unexcited_envelope_never_passes(self):
        tracker = StopTracker(source_end=0)
        for time in range(0, 102, 2):
            tracker.sample(time, {"slot": 0}, 0, R=0, T=0)
        self.assertFalse(tracker.converged)

    def test_zero_nodal_probe_allowed_when_envelope_excited(self):
        tracker = StopTracker(source_end=10)
        tracker.sample(0, {"node": 0}, 1, R=0.3, T=0.1)
        for time in range(12, 72, 2):
            tracker.sample(time, {"node": 0}, 1e-8, R=0.3, T=0.1)
        self.assertTrue(tracker.converged)

    def test_custom_decay_threshold_is_used(self):
        tracker = self.start(decay_threshold=1e-6)
        self.sample_tail(tracker)
        self.assertFalse(tracker.converged)

    def test_late_burst_revokes_convergence_and_requires_fresh_windows(self):
        tracker = self.sample_tail(self.start())
        self.assertTrue(tracker.converged)
        tracker.sample(72, {"slot": 0.2, "air": 0.1}, 1, R=0.3, T=0.1)
        self.assertFalse(tracker.converged)
        tracker.sample(74, {"slot": 1e-8, "air": 1e-8}, 1e-8, R=0.3, T=0.1)
        self.assertFalse(tracker.converged)

    def test_reject_invalid_values_and_changing_sample_contract(self):
        tracker = self.start()
        for time, probes, envelope in ((10, {"slot": 0, "air": 0}, 0),
                                      (12, {"slot": -1, "air": 0}, 0),
                                      (12, {"slot": 0}, 0),
                                      (12, {"slot": 0, "air": 0}, np.nan)):
            with self.assertRaises(ValueError):
                tracker.sample(time, probes, envelope, R=0.3, T=0.1)
        with self.assertRaises(ValueError):
            tracker.sample(12, {"slot": 0, "air": 0}, 0, source_end=11)


class DirectionMetricsTests(unittest.TestCase):
    def test_no_errors_means_unverified_even_large_contrast(self):
        result = direction_metrics(0.8, 0.2)
        self.assertFalse(result["verified"])
        self.assertIsNone(result["D"])
        self.assertAlmostEqual(result["raw_D"], 0.6)

    def test_five_error_floor_prevents_dividing_by_small_signal(self):
        result = direction_metrics(0.001, 0.0001, 0.0002, 0.0001)
        self.assertEqual(result["denominator_floor"], 0.0015)
        self.assertFalse(result["verified"])
        self.assertIsNone(result["D"])

    def test_finite_error_bound_and_angle_reversal(self):
        result = direction_metrics(0.8, 0.2, 0.01, 0.02)
        reverse = direction_metrics(0.2, 0.8, 0.02, 0.01)
        self.assertTrue(result["verified"])
        self.assertAlmostEqual(result["D"], 0.6)
        self.assertEqual(result["D"], -reverse["D"])
        self.assertEqual(result["error_D"], reverse["error_D"])
        for da in (-0.01, 0.01):
            for db in (-0.02, 0.02):
                perturbed = (0.8 + da - 0.2 - db) / (1 + da + db)
                self.assertLessEqual(abs(perturbed - result["D"]), result["error_D"])

    def test_coarse_grid_width_is_explicitly_unresolved(self):
        result = direction_metrics(0.8, 0.2, 0.01, 0.01, coarse=True)
        self.assertIsNone(result["FWHM"])
        self.assertTrue(result["unresolved_FWHM"])
        self.assertEqual(result["FWHM_status"], "unresolved_coarse_grid")
        json.dumps(result, allow_nan=False)

    def test_zero_signal_does_not_pass_even_zero_error(self):
        result = direction_metrics(0, 0, 0, 0)
        self.assertFalse(result["verified"])
        self.assertIsNone(result["raw_D"])

    def test_invalid_inputs(self):
        for args in ((-1, 0, 0, 0), (1, 0, -1, 0), (np.nan, 1, 0, 0)):
            with self.assertRaises(ValueError):
                direction_metrics(*args)


if __name__ == "__main__":
    unittest.main()
