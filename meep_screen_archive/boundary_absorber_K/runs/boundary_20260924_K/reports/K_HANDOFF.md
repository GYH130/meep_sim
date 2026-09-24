# K absorber verification — two cases only

This is not the 20-condition size study. No production gate is unlocked.

Only outer boundary implementation differs from PML12; frozen material and scientific tolerances are unchanged.

| Condition | Status |
|---|---|
| K_flat32p | NUMERICAL_FAILED |
| K_slot32p | NOT_RUN |

Full observables and fixed acceptance checks: per-task result.json.
Missing production evidence: mesh, stricter stop, mirror/straight-slot, s polarization and L48 checks.
Solver cutoff: 2026-09-24 20:30 China time. Rental expires approximately 21:53.
4 MPI ranks, serial only. K cap 10 h within original 48 h; no automatic retries.
Raw fields: separate checksum bundle, not GitHub. AI/RMB cost: unavailable, not inferred.
