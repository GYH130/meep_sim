# K — finite boundary verification, not a size sweep

Purpose: test Meep Absorber, 12 um, against the diagnosed PML12 setup without
changing the frozen Ordal material, geometry, source, r32 grid or numerical gates.
The 38% reduction seen in J is residual field intensity, NOT an emissivity error.

- New matched empty reference + flat p30: cap 4 h. Baseline slot p30: cap 6 h,
  launched only if the flat case and reference pass every original single-case gate.
- Serial 4 MPI ranks, single thread/rank, 70% measured container-available memory.
- K solver union cap 10 h within the original cumulative 48 h; no retries.
- User-confirmed rental expiry: approximately 2026-09-24 21:53 China time.
  Solver hard stop: 20:30; 83 min retained for preservation. Remaining budget may shorten any task.
- New boundary-aware reference identity; no PML reference reuse. Fresh fields, not timestep resume.
- Full target-frequency R/T/A, independent absorption, orders and stop evidence.
  A passing single case does not certify mesh/mirror/s-p/strict-stop or the 20-condition study.
- All data/code/cache stay under this data-disk directory; old files remain untouched.
- Background queue and backup scripts make no AI calls. RMB/token usage is unknown, not estimated.

Server commands (existing environment only):

```sh
cd /root/autodl-tmp/meep_sim/meep_screen/boundary_absorber_K
/root/miniconda3/envs/meep_sim/bin/python -B run_k.py prepare
/root/miniconda3/envs/meep_sim/bin/python -B run_k.py run
```

Status: `runs/boundary_20260924_K/status.json`; reports in that run's `reports/`.
No automatic GitHub push from the server. The separately launched local collector
copies raw evidence with checksums and updates only the authorized archive branch.
The Mac must remain awake and connected for collection; solver execution is independent.

Reference: https://meep.readthedocs.io/en/latest/FAQ/#why-are-the-fields-not-being-absorbed-by-the-pml
