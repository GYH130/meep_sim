"""Pure offline ownership/readback regressions; no Meep, MPI, or FDTD runs."""

import copy
import json
import unittest
from contextlib import nullcontext
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from ti2d.flux_reference import (
    assert_loaded_reference, buffer_check, layout_signature, make_partition, partition_spec,
)


class PartitionTests(unittest.TestCase):
    def test_explicit_owner_tree_is_json_portable_and_converted_to_meep_axes(self):
        expected = [["Y", 0.0], [["X", 0.0], 0, 1], [["X", 0.0], 2, 3]]
        self.assertEqual(json.loads(json.dumps(partition_spec(4))), expected)
        mp = SimpleNamespace(X=10, Y=20, BinaryPartition=lambda **kwargs: kwargs)
        self.assertEqual(make_partition(mp, 4),
                         {"data": [(20, 0.0), [(10, 0.0), 0, 1], [(10, 0.0), 2, 3]]})
        spec = partition_spec(4)
        spec[1][1] = 2
        self.assertEqual(partition_spec(4), expected)

    def test_only_exact_four_integer_ranks_are_accepted(self):
        for ranks in (1, 2, 8, 4.0, "4", True, np.bool_(True), None):
            with self.subTest(ranks=ranks), self.assertRaises(ValueError):
                partition_spec(ranks)
        self.assertEqual(partition_spec(np.int64(4)), partition_spec(4))


class LayoutSignatureTests(unittest.TestCase):
    def datasets(self):
        return {"gv_nums": np.array([10, 30, 0, 10, 30, 0, 10, 30, 0, 10, 30, 0]),
                "gv_origins": np.array([-1., -3., 0., 0., -3., 0., -1., 0., 0., 0., 0., 0.])}

    def signature(self, datasets=None, owners=(0, 1, 2, 3)):
        datasets = self.datasets() if datasets is None else datasets
        fake_h5py = SimpleNamespace(File=lambda filename, mode: nullcontext(datasets))
        with patch.dict("sys.modules", {"h5py": fake_h5py}):
            return layout_signature("fake-chunks.h5", owners)

    def test_owners_are_part_of_exact_signature(self):
        before = self.signature()
        after = self.signature(owners=(2, 3, 0, 1))
        self.assertEqual(before["datasets"], after["datasets"])
        self.assertNotEqual(before["sha256"], after["sha256"])
        self.assertEqual(before["owners"], [0, 1, 2, 3])
        json.dumps(before, allow_nan=False)

    def test_logical_values_ignore_hdf5_storage_integer_width_and_byte_order(self):
        before = self.signature()
        data = self.datasets()
        data["gv_nums"] = data["gv_nums"].astype(">i4")
        data["gv_origins"] = data["gv_origins"].astype(">f8")
        self.assertEqual(self.signature(data), before)
        data["gv_nums"][0] += 1
        self.assertNotEqual(self.signature(data)["sha256"], before["sha256"])

    def test_coordinate_change_is_detected(self):
        before = self.signature()
        data = self.datasets()
        data["gv_origins"][0] += 1e-10
        self.assertNotEqual(self.signature(data)["sha256"], before["sha256"])

    def test_invalid_or_missing_layout_and_owner_information_is_rejected(self):
        for owners in ([], [0., 1., 2., 3.], [0, 1, -1, 3], [[0, 1], [2, 3]], [True] * 4):
            with self.subTest(owners=owners), self.assertRaises(ValueError):
                self.signature(owners=owners)
        for change in ("missing", "size", "nonfinite", "noninteger", "negative"):
            data = self.datasets()
            if change == "missing":
                del data["gv_nums"]
            elif change == "size":
                data["gv_nums"] = np.zeros(9)
            elif change == "nonfinite":
                data["gv_origins"][0] = np.nan
            elif change == "noninteger":
                data["gv_nums"] = np.full(12, 1.5)
            elif change == "negative":
                data["gv_nums"][0] = -1
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.signature(data)


class BufferCheckTests(unittest.TestCase):
    def pair(self, value=1):
        return {"E": np.array([value, 2j * value]), "H": np.array([3 * value, -1j * value])}

    def test_matching_complex_buffers_and_empty_owner_are_valid(self):
        wanted = self.pair()
        check = buffer_check(SimpleNamespace(**wanted), (wanted["E"].copy(), wanted["H"].copy()))
        self.assertTrue(check["valid"])
        self.assertTrue(check["shape_match"])
        self.assertTrue(check["finite"])
        self.assertEqual(check["expected_nonzero"], 4)
        self.assertEqual(check["loaded_power"], 15)
        for zero in (self.pair(0), {"E": np.array([]), "H": np.array([])}):
            check = buffer_check(zero, zero)
            self.assertTrue(check["valid"])
            self.assertEqual(check["loaded_power"], 0)
        json.dumps(check, allow_nan=False)

    def test_shape_nonfinite_and_wrong_sign_fail_with_reports(self):
        wanted = self.pair()
        for change in ("shape", "nan", "wrong_sign", "missing", "text"):
            actual = self.pair()
            if change == "shape":
                actual["E"] = actual["E"][:, None]
            elif change == "nan":
                actual["H"][0] = np.nan
            elif change == "wrong_sign":
                actual["E"] *= -1
            elif change == "missing":
                del actual["H"]
            elif change == "text":
                actual["H"] = np.array(["a", "b"])
            with self.subTest(change=change):
                check = buffer_check(wanted, actual)
                self.assertFalse(check["valid"])
                self.assertTrue(check["errors"])
                json.dumps(check, allow_nan=False)

    def test_tolerances_scale_with_the_actual_signal(self):
        for scale in (1e-100, 1, 1e100):
            wanted = self.pair(scale)
            actual = {name: value * (1 + 1e-13) for name, value in wanted.items()}
            self.assertTrue(buffer_check(wanted, actual)["valid"])
            actual = {name: np.zeros_like(value) for name, value in wanted.items()}
            self.assertFalse(buffer_check(wanted, actual)["valid"])
            actual = {name: value * (1 + 1e-8) for name, value in wanted.items()}
            self.assertFalse(buffer_check(wanted, actual)["valid"])

    def test_large_integer_power_does_not_wrap_and_overflow_is_reported(self):
        ints = {"E": np.array([2 ** 62], dtype=np.int64), "H": np.array([0], dtype=np.int64)}
        check = buffer_check(ints, ints)
        self.assertTrue(check["valid"])
        self.assertEqual(check["loaded_power"], float(2 ** 124))
        huge = self.pair(1e308)
        check = buffer_check(huge, huge)
        self.assertFalse(check["valid"])
        json.dumps(check, allow_nan=False)


class CollectiveReadbackTests(unittest.TestCase):
    def buffers(self):
        result = []
        for rank in range(4):
            E = np.zeros(4, dtype=complex)
            H = np.zeros(4, dtype=complex)
            if rank < 2:
                E[rank * 2:rank * 2 + 2] = 1 + rank + 1j
                H[rank * 2:rank * 2 + 2] = 2 + rank - 1j
            result.append({"E": E, "H": H})
        return result

    def check_all_ranks(self, expected, loaded, *, power=2.0, flux=-2.0, error=None):
        reports = []
        for rank in range(4):
            check = buffer_check({name: -value for name, value in expected[rank].items()}, loaded[rank])
            reports.append({"rank": rank, "input_power": power, "normalized_flux": None,
                            "buffer_check": check,
                            "errors": [] if check["valid"] else ["REFERENCE_LOADED_BUFFER_MISMATCH"]})
        audits = []
        for rank in range(4):
            calls = []
            def readback(monitor):
                calls.append("data")
                return SimpleNamespace(**loaded[rank])
            def fluxes(monitor):
                calls.append("flux")
                return [flux]
            def gather(local):
                calls.append("gather")
                gathered = copy.deepcopy(reports)
                gathered[rank] = local
                return gathered
            sim = SimpleNamespace(get_flux_data=readback)
            mp = SimpleNamespace(get_fluxes=fluxes)
            comm = SimpleNamespace(rank=rank, size=4, allgather=gather)
            if error:
                with self.assertRaisesRegex(ValueError, error) as caught:
                    assert_loaded_reference(sim, object(), expected[rank], power, mp, comm)
                audit = caught.exception.audit
                self.assertFalse(audit["passed"])
            else:
                audit = assert_loaded_reference(sim, object(), expected[rank], power, mp, comm)
                self.assertTrue(audit["passed"])
            self.assertEqual(calls, ["data", "flux", "gather"])
            json.dumps(audit, allow_nan=False)
            audits.append(audit)
        return audits

    def test_fixed_partition_preserves_nonzero_rank_buffers_and_passes(self):
        expected = self.buffers()
        loaded = [{name: -value for name, value in pair.items()} for pair in expected]
        audits = self.check_all_ranks(expected, loaded)
        self.assertGreater(audits[0]["global_buffer_power"]["loaded_power"], 0)
        self.assertEqual(audits[2]["ranks"][2]["buffer_check"]["loaded_power"], 0)
        self.assertEqual(audits[0]["ranks"][0]["normalized_flux"], -1)

    def test_old_rank_map_shift_loses_buffers_and_fails_on_every_rank(self):
        expected = self.buffers()
        loaded = [{name: np.zeros_like(value) for name, value in pair.items()} for pair in expected]
        # Cached active owners0/1 -> target owners2/3: local loads retain zeros.
        self.check_all_ranks(expected, loaded, flux=0, error="REFERENCE_LOADED_AUDIT_FAILED")

    def test_rank_permutation_fails_even_if_global_power_and_flux_look_correct(self):
        expected = self.buffers()
        loaded = [{name: -value for name, value in expected[(rank + 2) % 4].items()}
                  for rank in range(4)]
        self.check_all_ranks(expected, loaded, error="REFERENCE_LOADED_BUFFER_MISMATCH")

    def test_zero_reference_is_rejected_even_with_apparently_correct_flux(self):
        zeros = [{"E": np.zeros(4), "H": np.zeros(4)} for _ in range(4)]
        self.check_all_ranks(zeros, zeros, error="REFERENCE_GLOBAL_EXPECTED_POWER_NOT_POSITIVE")

    def test_wrong_flux_sign_and_magnitude_are_rejected(self):
        expected = self.buffers()
        loaded = [{name: -value for name, value in pair.items()} for pair in expected]
        for flux in (2., -1., -2. + 1e-7, np.nan):
            with self.subTest(flux=flux):
                self.check_all_ranks(expected, loaded, flux=flux, error="REFERENCE_LOADED_AUDIT_FAILED")

    def test_invalid_input_does_not_skip_collective_readbacks(self):
        expected = self.buffers()
        loaded = [{name: -value for name, value in pair.items()} for pair in expected]
        for power in (0, -1, np.nan, np.inf, True, 2j, None):
            # All fake remote reports must carry serializable invalid markers,
            # matching the actual helper's validated per-rank representation.
            calls = []
            sim = SimpleNamespace(get_flux_data=lambda monitor: (calls.append("data") or loaded[0]))
            mp = SimpleNamespace(get_fluxes=lambda monitor: (calls.append("flux") or [-2.]))
            def gather(local):
                calls.append("gather")
                return [dict(copy.deepcopy(local), rank=rank) for rank in range(4)]
            comm = SimpleNamespace(rank=0, size=4, allgather=gather)
            with self.subTest(power=power), self.assertRaisesRegex(ValueError, "REFERENCE_INPUT_POWER_INVALID"):
                assert_loaded_reference(sim, object(), expected[0], power, mp, comm)
            self.assertEqual(calls, ["data", "flux", "gather"])

    def test_local_readback_exception_is_gathered_before_collective_failure(self):
        calls = []
        def readback(monitor):
            calls.append("data")
            raise AttributeError("simulated rank-local readback error")
        def gather(local):
            calls.append("gather")
            return [dict(copy.deepcopy(local), rank=rank) for rank in range(4)]
        sim = SimpleNamespace(get_flux_data=readback)
        mp = SimpleNamespace(get_fluxes=lambda monitor: (calls.append("flux") or [-2.]))
        comm = SimpleNamespace(rank=0, size=4, allgather=gather)
        with self.assertRaisesRegex(ValueError, "REFERENCE_BUFFER_READBACK_FAILED"):
            assert_loaded_reference(sim, object(), self.buffers()[0], 2., mp, comm)
        self.assertEqual(calls, ["data", "flux", "gather"])


if __name__ == "__main__":
    unittest.main()
