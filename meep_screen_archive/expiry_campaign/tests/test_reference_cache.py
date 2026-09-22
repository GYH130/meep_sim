"""Cache corruption and path validation using temporary bytes, never Meep."""
from copy import deepcopy
from pathlib import Path
import tempfile
import unittest

from ti2d.common import digest_file, digest_object
from ti2d.reference_cache import validated_reference


class ReferenceCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.attempt = "1" * 32
        directory = self.root / self.attempt
        directory.mkdir()
        names = ["chunks.h5", "reference_planes.npz", *(f"flux_rank{rank:03d}.npz" for rank in range(4))]
        evidence = {}
        for name in names:
            path = directory / name
            path.write_bytes(("synthetic evidence " + name).encode())
            evidence[f"{self.attempt}/{name}"] = digest_file(path)
        self.identity = {"mpi_ranks": 4, "implementation_sha256": "current-code", "fixed_partition": "fixed"}
        signature = {"schema_version": 1, "owners": [0, 1, 2, 3], "datasets": {
            "gv_nums": {"shape": [12], "values": [10, 20, 0] * 4},
            "gv_origins": {"shape": [12], "values": [-1., -2., 0.] * 4}}}
        signature["sha256"] = digest_object(signature)
        self.candidate = {"status": "QUALIFIED", "identity": deepcopy(self.identity),
            "checks": {name: True for name in ("stop", "downward_flux", "direction", "plane_closure",
                                                "empty_pseudo_reflection", "empty_transmission")},
            "stop": {"converged": True, "stop_reason": "CONVERGED", "source_on_non_pml_envelope_peak": .02},
            "input_power": .5, "pseudo_reflection": .00001, "chunk_layout_signature": signature,
            "attempt": self.attempt, "evidence_sha256": evidence}

    def valid(self, candidate=None):
        return validated_reference(self.candidate if candidate is None else candidate,
                                   self.root, self.identity, 4)

    def test_complete_matching_candidate_passes_without_mutation(self):
        before = deepcopy(self.candidate)
        self.assertTrue(self.valid())
        self.assertEqual(self.candidate, before)

    def test_identity_qualification_and_stop_must_match(self):
        changes = [("status", "NUMERICAL_FAILED"), ("identity", {}), ("checks", {}),
                   ("checks", {**self.candidate["checks"], "stop": 1}), ("stop", None),
                   ("stop", {**self.candidate["stop"], "stop_reason": "TIME_LIMIT_UNQUALIFIED"}),
                   ("schema_version", True)]
        for key, value in changes:
            candidate = deepcopy(self.candidate)
            candidate[key] = value
            with self.subTest(key=key, value=value):
                self.assertFalse(self.valid(candidate))

    def test_reference_anchor_and_power_are_positive_finite_real_scalars(self):
        for value in (None, 0, -1, True, "1", float("nan"), float("inf"), 1j):
            for field in ("input_power", "source_on_non_pml_envelope_peak"):
                candidate = deepcopy(self.candidate)
                destination = candidate["stop"] if field.startswith("source_on") else candidate
                destination[field] = value
                with self.subTest(field=field, value=value):
                    self.assertFalse(self.valid(candidate))

    def test_empty_partial_extra_or_wrong_hash_evidence_is_rejected(self):
        for change in ("empty", "missing_rank", "extra", "wrong_hash", "bad_type"):
            candidate = deepcopy(self.candidate)
            if change == "empty": candidate["evidence_sha256"] = {}
            elif change == "missing_rank": candidate["evidence_sha256"].pop(f"{self.attempt}/flux_rank003.npz")
            elif change == "extra": candidate["evidence_sha256"]["unrelated.txt"] = "0" * 64
            elif change == "wrong_hash": candidate["evidence_sha256"][f"{self.attempt}/chunks.h5"] = "0" * 64
            else: candidate["evidence_sha256"] = []
            with self.subTest(change=change):
                self.assertFalse(self.valid(candidate))

    def test_missing_or_changed_actual_file_is_rejected(self):
        path = self.root / self.attempt / "flux_rank003.npz"
        path.write_bytes(b"damaged")
        self.assertFalse(self.valid())
        path.unlink()
        self.assertFalse(self.valid())

    def test_attempt_and_evidence_paths_cannot_escape(self):
        for attempt in ("../" + self.attempt, "/" + self.attempt, ".", "", "x" * 32, self.attempt + "/child"):
            candidate = deepcopy(self.candidate)
            candidate["attempt"] = attempt
            with self.subTest(attempt=attempt):
                self.assertFalse(self.valid(candidate))
        for path in ("/tmp/flux_rank000.npz", f"{self.attempt}/../flux_rank000.npz"):
            candidate = deepcopy(self.candidate)
            value = candidate["evidence_sha256"].pop(f"{self.attempt}/flux_rank000.npz")
            candidate["evidence_sha256"][path] = value
            self.assertFalse(self.valid(candidate))

    def test_file_symlink_is_rejected_even_with_correct_target_hash(self):
        original = self.root / self.attempt / "flux_rank000.npz"
        external = self.root / "external.npz"
        original.rename(external)
        original.symlink_to(external)
        self.assertFalse(self.valid())

    def test_layout_structure_hash_and_owners_are_checked(self):
        for change in ("missing", "wrong_owner", "shape", "nonfinite", "stale_hash", "bad_counts"):
            candidate = deepcopy(self.candidate)
            signature = candidate["chunk_layout_signature"]
            if change == "missing": signature.pop("datasets")
            elif change == "wrong_owner": signature["owners"][0] = 4
            elif change == "shape": signature["datasets"]["gv_nums"]["shape"] = [11]
            elif change == "nonfinite": signature["datasets"]["gv_origins"]["values"][0] = float("nan")
            elif change == "stale_hash": signature["datasets"]["gv_origins"]["values"][0] += .25
            else: signature["datasets"]["gv_nums"]["values"][0] = 1.5
            with self.subTest(change=change):
                self.assertFalse(self.valid(candidate))

    def test_malformed_candidates_return_false_without_raising(self):
        for candidate in (None, [], True, 42, "text", {"status": "QUALIFIED", "identity": None}):
            self.assertFalse(validated_reference(candidate, self.root, self.identity, 4))
        for ranks in (0, -1, True, 4.0, "4"):
            self.assertFalse(validated_reference(self.candidate, self.root, self.identity, ranks))


if __name__ == "__main__":
    unittest.main()
