# L: bounded homogeneous-boundary screen

This stage does not certify the slotted solver or start the 20-condition size study.

K stopped at its empty-reference gate: backward/forward power was
0.0017166551733512714, above the unchanged 0.001 limit. Its decay and direction
checks passed. Independent offline projection of saved K complex fields found
the same uniform backward wave at both planes, not a higher-order projection
artifact. The full-width planar metal and slot were not run in K.

L uses the unchanged K solver implementation and frozen Ordal material files.
Only the computational width of homogeneous air/slab cases is reduced to 2 um;
the real 50 um patterned geometry is not changed. The dummy cavity metadata is
never instantiated in these flat cases. First reproduce K's failed reference
ratio and incident power per unit width. A failed replay stops this queue.

After a successful diagnostic replay, test Absorber R_asymptotic=1e-8, then
1e-6 only if necessary. Each candidate retains all reference and flat-slab
acceptance checks; selection also requires reference reflection below 0.0005.
No scaling or clipping of measured powers is permitted. A selected candidate
is not evidence for grazing slot diffraction orders or full-cell qualification.

Maximum three solves, 20 minutes each, one hour total solver union, four MPI
ranks, one task at a time, no retries. All time counts in the shared 48-hour
ledger. Hard solver cutoff remains 2026-09-24 20:30 China time. Full-width
configurations are recorded under configs/future, explicitly not queued because
the remaining rental window is insufficient. Stop after screening.

Entry points on the server, using the existing meep_sim environment:

- `python -B prepare_l.py`: single-use input/backup/fingerprint preparation.
- `python -B -m unittest discover -s tests -q`: offline tests, no FDTD.
- `python -B -m unittest test_l_control -q`: offline controller tests.
- `python -B run_l.py`: single-use bounded queue; no automatic resume/retry.
- Read `runs/boundary_20260924_L/status.json` for a compact live status.
- Read `runs/boundary_20260924_L/reports/L_HANDOFF.md` after a checkpoint.
- `python -B export_l.py --label prepared`: export before launch only; do not
  export while a solver is writing evidence. The queue exports between tasks.

Immutable raw evidence is downloaded to the Mac. Sanitized compact snapshots
go only to the previously authorized archive branch, never to main. Large
arrays, cache and full logs are not included in the public compact archive.
The local collector requires this Mac to stay awake and connected. The queue
itself makes no AI/API calls. Actual AI billing is unavailable, not zero.
