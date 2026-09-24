"""Pure validation of a complete, qualified incident-field cache candidate.

No Meep or MPI calls. Invalid, incomplete, unsafe, and unreadable candidates
return False so the caller can collectively choose a fresh reference solve.
"""
from __future__ import annotations

import math
from pathlib import Path
import re

from .common import digest_file, digest_object


_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_ATTEMPT = re.compile(r"[0-9a-f]{32}\Z")
_CHECKS = {"stop", "downward_flux", "direction", "plane_closure",
           "empty_pseudo_reflection", "empty_transmission"}


def _number(value, *, positive=False):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and (not positive or value > 0))


def _valid_layout(signature, ranks):
    if not isinstance(signature, dict) or set(signature) != {"schema_version", "datasets", "owners", "sha256"}:
        return False
    if type(signature["schema_version"]) is not int or signature["schema_version"] != 1:
        return False
    owners = signature["owners"]
    if (not isinstance(owners, list) or not owners or
            any(type(owner) is not int or not 0 <= owner < ranks for owner in owners)):
        return False
    datasets = signature["datasets"]
    if not isinstance(datasets, dict) or set(datasets) != {"gv_nums", "gv_origins"}:
        return False
    count = 3 * len(owners)
    for name, dataset in datasets.items():
        if not isinstance(dataset, dict) or set(dataset) != {"shape", "values"}:
            return False
        shape, values = dataset["shape"], dataset["values"]
        if (not isinstance(shape, list) or len(shape) != 1 or type(shape[0]) is not int or
                shape[0] != count or not isinstance(values, list) or len(values) != count):
            return False
        if name == "gv_nums":
            if any(type(value) is not int or value < 0 for value in values):
                return False
        elif not all(_number(value) for value in values):
            return False
    expected = digest_object({key: value for key, value in signature.items() if key != "sha256"})
    return isinstance(signature["sha256"], str) and signature["sha256"] == expected


def validated_reference(candidate, cache_dir, expected_identity, ranks):
    """Return whether all mandatory cache evidence and metadata are trustworthy.

    Each UUID attempt must contain exactly the six hashed data files for four
    ranks (generalized to ``ranks + 2``). Paths must be direct regular files in
    that attempt; neither traversal, external paths, nor symlinks are accepted.
    The manifest itself is intentionally not part of its own evidence hash map.
    """
    try:
        if type(ranks) is not int or ranks < 1:
            return False
        if not isinstance(candidate, dict) or not isinstance(expected_identity, dict):
            return False
        if (candidate.get("status") != "QUALIFIED" or not isinstance(candidate.get("identity"), dict)
                or candidate["identity"] != expected_identity or expected_identity.get("mpi_ranks") != ranks):
            return False
        if "schema_version" in candidate and (type(candidate["schema_version"]) is not int or
                                               candidate["schema_version"] != 1):
            return False
        checks = candidate.get("checks")
        if not isinstance(checks, dict) or not _CHECKS.issubset(checks) or not all(value is True for value in checks.values()):
            return False
        stop = candidate.get("stop")
        if (not isinstance(stop, dict) or stop.get("converged") is not True
                or stop.get("stop_reason") != "CONVERGED"
                or not _number(stop.get("source_on_non_pml_envelope_peak"), positive=True)
                or not _number(candidate.get("input_power"), positive=True)):
            return False
        pseudo = candidate.get("pseudo_reflection")
        if not _number(pseudo) or not 0 <= pseudo < .001:
            return False
        if not _valid_layout(candidate.get("chunk_layout_signature"), ranks):
            return False
        attempt = candidate.get("attempt")
        if not isinstance(attempt, str) or not _ATTEMPT.fullmatch(attempt):
            return False
        root = Path(cache_dir).resolve(strict=True)
        attempt_dir = root / attempt
        if attempt_dir.is_symlink() or not attempt_dir.is_dir() or attempt_dir.resolve(strict=True).parent != root:
            return False
        names = {"chunks.h5", "reference_planes.npz", *(f"flux_rank{rank:03d}.npz" for rank in range(ranks))}
        required = {f"{attempt}/{name}" for name in names}
        evidence = candidate.get("evidence_sha256")
        if not isinstance(evidence, dict) or set(evidence) != required:
            return False
        for relative, expected_hash in evidence.items():
            if not isinstance(expected_hash, str) or not _SHA256.fullmatch(expected_hash):
                return False
            path = root / relative
            if path.is_symlink() or not path.is_file() or path.resolve(strict=True).parent != attempt_dir:
                return False
            if digest_file(path) != expected_hash:
                return False
        return True
    except (OSError, ValueError, TypeError, KeyError, AttributeError, OverflowError, RuntimeError):
        return False
