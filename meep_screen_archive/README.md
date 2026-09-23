# Ti 2D stage archive — incomplete study

Pure Ti, specified frozen Ordal Route A optical model; 10.5 µm, 2D only.
This is NOT a validated size–performance study, 3D-hole prediction or real high-temperature sample measurement.

## Evidence classification
- `VERIFICATION.json` lists individually qualified cases with checked raw-evidence SHA256.
- Other results remain failures, time-limited or diagnostic-only; no automatic upgrade.
- Old flat r32 with smoothing passed; old slanted r32 overflowed.
- Smoothing-off short test reached t≈98.83 before its wall limit; source ends at 525. It does NOT prove full stability.
- H is a 10-condition validation subset with smoothing off, serial 4 MPI ranks, fixed gates, no retries or size scan.
- A failed required condition stops expansion. Missing checks stay missing. Original per-solver cap remains 6h. The former September 23 16:52 China rental deadline was withdrawn by the user; exact September 24 expiry time is not confirmed.

## Latest bounded tail diagnostics
- I stopped at its 2.5h cap (t914.625); original status and initial report are preserved for audit, not a completed study.
- I `reports_recovered/` is the corrected offline analysis. A single terminal partial interval is explicitly excluded from FFT; original arrays remain unchanged.
- J is a NEW zero-field matched pair: PML6 then PML12, target t1200, individual caps 4h and 5h, shared stage cap 9h within cumulative 48h.
- PML12 requires completed control and the unchanged replay-through-t1000 gate. No automatic retries or production scan.
- J status is in `expiry_campaign/runs/tail_20260923_J/status.json`; reports, when available, in its `reports_recovered/` directory. PREPARED/RUNNING is not completed.
- Tail spectra and field amplitudes do not constitute qualified optical reflectance, emission or size trends.

## Preservation and costs
This compact archive contains code, complete material inputs, configurations, small results/reports and raw-field fingerprints.
Large NPZ/HDF5 fields, caches, credentials, full conversation and raw logs are deliberately NOT committed.
Full pre-expiry raw evidence was separately downloaded as `preexpiry_evidence_20260922.tar.gz` (SHA256 17b1bc4eae4388bdcd1ccfe15496097ce20305102d0d3e5de1ab4eb8b1a24573).
An independent post-run raw bundle is generated on the server. Do not delete/release the instance before retrieving it.
`budget_ledger.json` measures union solver time and ranks; AI usage/billing and RMB costs are unavailable, not inferred. Background queue has no AI/API calls.

## Running / recovery
Remote stage `/root/autodl-tmp/meep_sim/meep_screen`, existing `meep_sim` Conda only.
See `expiry_campaign/HANDOFF.md`. Preserve absolute inputs or regenerate a NEW hashed manifest for another server; do not modify an already locked manifest.
This archive is a snapshot in a NEW GitHub branch. Main and the existing research branch are unchanged.
