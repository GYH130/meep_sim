"""Offline intensity-trace regressions; no Meep import or FDTD execution."""

import json
import unittest

import numpy as np

from ti2d.diagnostics import StopTracker


class ReferenceAnchoredStopTests(unittest.TestCase):
    reference_peak = 0.976
    back_peak = 1.98e-9
    back_tail = 0.97 * back_peak
    envelope_tail = 2.6e-9 * reference_peak

    def tracker(self, **kwargs):
        options = dict(source_end=10, reference_intensity_peak=self.reference_peak,
                       reference_provenance="synthetic reference source-on |E|^2 peak; not fitted to tail")
        options.update(kwargs)
        return StopTracker(**options)

    def excite(self, tracker, scale=1):
        tracker.sample(0, {"front_air": 0.4 * scale, "back_air": self.back_peak * scale},
                       self.reference_peak * scale, R=0.5, T=1e-9)
        tracker.sample(10, {"front_air": 0.1 * scale, "back_air": self.back_peak * scale},
                       0.2 * scale, R=0.5, T=1e-9)

    def tail(self, tracker, *, end=70, back=None, envelope=None, R=0.5, T=1e-9, scale=1):
        back = self.back_tail if back is None else back
        envelope = self.envelope_tail if envelope is None else envelope
        for time in range(12, end + 1, 2):
            tracker.sample(time, {"front_air": 1e-9 * scale, "back_air": back * scale},
                           envelope * scale, R=R, T=T)
        return tracker

    def test_saved_trace_scale_replay_dark_back_probe_passes_both_thresholds(self):
        # Reproduces the diagnosed scales, not a claim to contain the saved run:
        # back historical peak ~1.98e-9, local tail/peak .97, envelope peak .976,
        # non-PML tail/peak 2.6e-9. The legacy predicate rejects every tail sample.
        self.assertGreater(self.back_tail / self.back_peak, 1e-4)
        for threshold in (1e-4, 1e-6):
            with self.subTest(decay_threshold=threshold):
                tracker = self.tracker(decay_threshold=threshold)
                self.excite(tracker)
                self.tail(tracker, end=68)
                self.assertFalse(tracker.converged)
                tracker.sample(70, {"front_air": 1e-9, "back_air": self.back_tail},
                               self.envelope_tail, R=0.5, T=1e-9)
                self.assertTrue(tracker.converged)
                sample = tracker.trace[-1]
                self.assertEqual(sample["probe_peak"]["back_air"], self.back_peak)
                self.assertEqual(sample["probe_normalization_denominator"]["back_air"], self.reference_peak)
                self.assertAlmostEqual(sample["probe_ratio"]["back_air"], self.back_tail / self.reference_peak)
                self.assertEqual(sample["probe_classification"]["back_air"], "below_reference_decay_budget")
                self.assertEqual(sample["probe_normalization_source"]["back_air"], "reference_excitation_floor")
                self.assertTrue(sample["normalization"]["all_probes_required"])
                self.assertEqual(tracker.summary()["required_windows"], 3)
                self.assertEqual(tracker.summary()["window"], 20)
                self.assertEqual(tracker.summary()["rt_tolerance"], 1e-3)
                json.dumps(tracker.trace, allow_nan=False)
                json.dumps(tracker.summary(), allow_nan=False)

    def test_meaningful_persistent_dark_probe_residual_is_not_dropped(self):
        for threshold in (1e-4, 1e-6):
            tracker = self.tracker(decay_threshold=threshold)
            self.excite(tracker)
            self.tail(tracker, end=150, back=2 * threshold * self.reference_peak)
            self.assertFalse(tracker.converged)
            self.assertGreater(tracker.trace[-1]["probe_ratio"]["back_air"], threshold)

    def test_sustained_spurious_high_field_is_not_hidden_by_floor(self):
        tracker = self.tracker()
        self.excite(tracker)
        self.tail(tracker, end=150, back=100 * self.reference_peak,
                  envelope=100 * self.reference_peak)
        self.assertFalse(tracker.converged)
        self.assertEqual(tracker.trace[-1]["probe_ratio"]["back_air"], 1)
        self.assertEqual(tracker.trace[-1]["envelope_ratio"], 1)
        self.assertEqual(tracker.summary()["normalization"]["reference_intensity_peak"], self.reference_peak)

    def test_strict_threshold_still_rejects_residual_between_thresholds(self):
        ordinary, strict = self.tracker(), self.tracker(decay_threshold=1e-6)
        for tracker in (ordinary, strict):
            self.excite(tracker)
            self.tail(tracker, back=1e-5 * self.reference_peak)
        self.assertTrue(ordinary.converged)
        self.assertFalse(strict.converged)

    def test_default_reference_run_anchors_to_source_on_envelope_and_freezes(self):
        tracker = StopTracker(source_end=10)
        self.excite(tracker)
        self.tail(tracker)
        self.assertTrue(tracker.converged)
        normalization = tracker.summary()["normalization"]
        self.assertEqual(normalization["reference_peak_source"], "source_on_non_pml_envelope")
        self.assertEqual(normalization["denominator_floor"], self.reference_peak)
        self.assertTrue(normalization["reference_frozen"])
        tracker.sample(72, {"front_air": 2, "back_air": 2}, 10, R=0.5, T=1e-9)
        self.assertFalse(tracker.converged)
        self.assertEqual(tracker.summary()["normalization"]["denominator_floor"], self.reference_peak)

    def test_tail_cannot_create_missing_reference_excitation(self):
        tracker = StopTracker(source_end=10)
        tracker.sample(0, {"front_air": 0, "back_air": 0}, 0, R=0.5, T=1e-9)
        self.tail(tracker, end=150)
        self.assertFalse(tracker.converged)
        self.assertFalse(tracker.summary()["normalization"]["reference_ready"])
        self.assertIsNone(tracker.trace[-1]["probe_ratio"]["back_air"])
        self.assertIsNone(tracker.trace[-1]["envelope_ratio"])

    def test_floor_is_intensity_scale_invariant(self):
        # Includes small physical intensities; no dimensionful magic epsilon.
        for scale in (1e-12, 1, 1e12):
            tracker = self.tracker(reference_intensity_peak=self.reference_peak * scale)
            self.excite(tracker, scale=scale)
            self.tail(tracker, scale=scale)
            self.assertTrue(tracker.converged)
            self.assertAlmostEqual(tracker.trace[-1]["probe_ratio"]["back_air"],
                                   self.back_tail / self.reference_peak)

    def test_fraction_can_tighten_but_not_relax_reference_budget(self):
        tracker = self.tracker(reference_floor_fraction=0.01)
        self.excite(tracker)
        self.tail(tracker, back=1e-5 * self.reference_peak)
        self.assertFalse(tracker.converged)
        normalization = tracker.summary()["normalization"]
        self.assertEqual(normalization["denominator_floor"], 0.01 * self.reference_peak)
        self.assertAlmostEqual(normalization["weak_probe_residual_bound"],
                               1e-4 * 0.01 * self.reference_peak)

    def test_missing_either_rt_component_still_prevents_convergence(self):
        for R, T in ((None, None), (None, 1e-9), (0.5, None)):
            tracker = self.tracker()
            self.excite(tracker)
            self.tail(tracker, end=150, R=R, T=T)
            self.assertFalse(tracker.converged)

    def test_reference_provenance_and_intensity_contract_validation(self):
        for peak in (0, -1, np.nan, np.inf, True, 1 + 0j, [1]):
            with self.subTest(peak=peak), self.assertRaises(ValueError):
                self.tracker(reference_intensity_peak=peak)
        for fraction in (0, -1, 1.01, np.nan, np.inf, True, 1j, [1]):
            with self.subTest(fraction=fraction), self.assertRaises(ValueError):
                self.tracker(reference_floor_fraction=fraction)
        for provenance in (None, "", "  ", 123):
            with self.subTest(provenance=provenance), self.assertRaises(ValueError):
                self.tracker(reference_provenance=provenance)
        with self.assertRaises(ValueError):
            StopTracker(reference_provenance="missing the peak")
        with self.assertRaises(ValueError):
            self.tracker(reference_intensity_peak=np.nextafter(0., 1.), reference_floor_fraction=0.01)
        for windows in (0, -1, 1.5, np.nan, np.inf, True, np.bool_(True), 1j, [3]):
            with self.subTest(windows=windows), self.assertRaises(ValueError):
                self.tracker(consecutive_windows=windows)


if __name__ == "__main__":
    unittest.main()
