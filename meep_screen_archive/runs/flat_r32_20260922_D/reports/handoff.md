# Bounded Ti 2D stage handoff

State: **VALIDATION_FAILED_OR_INCOMPLETE**. Pilot conditions qualified: **0/20**.
Qualification gate: **NOT PASSED**.

Scope: specified Ordal Route A Ti optical model, 2D slanted slots, 10.5 um, +/-30 deg, s/p only. No oxide layer; real-sample oxidation and high-temperature validity remain unknown.

## Resource use

- Active-solver interval union: 8.3670 h / 48 h.
- MPI rank-hours: 33.4682; calendar span: 17.1570 h.
- No AI calls are made by this queue. AI usage and RMB costs are unavailable unless an actual account bill is supplied; no price is invented.
- Concurrency: one locked two-job test uses required validation conditions, each limited to half the group cap of 70% container-available memory. More than 30% slowdown or incomparable cache state returns to serial; heterogeneous later cases conservatively stay serial. Probe status: NOT_RUN. Restarting an interrupted condition is not time-step checkpoint recovery.

## Interpretation and unfinished work

- Missing locked mesh comparison
- Missing locked mirror comparison
- Missing locked strict comparison

- Preserve unclipped R, T, A_flux, A_vol. Algebraic R+T+A_flux=1 is not independent conservation evidence.
- A numeric value in conditions.csv does not imply qualification; use its status and this gate.
- Fixed surface fill with changed period changes opening width; fixed opening with changed tilt changes normal width. Only local finite differences are reported.
- A missing validated error bound means contrast is unverified, not zero uncertainty. Effects within paired error bounds are unresolved at current precision.
- No formal peak angle, angular spectrum, FWHM, optimal design, 3D or high-temperature conclusion.

## Evidence and next action

- `conditions.csv`, `paired_emissivity.csv`, `local_sensitivity.csv`, `numerical_error_evidence.csv`, `summary.json`, `figures/`. Observed grid/stop spreads are evidence, not a rigorous universal error bound.
- Original arrays and stop curves remain in task attempt directories; previous files are not overwritten.
- Inspect a failed condition's result.json and short log tail before authorizing any targeted repair or rerun. No automatic scan expansion.
- For status use `python -B -m ti2d.cli status --manifest /root/autodl-tmp/meep_sim/meep_screen/runs/flat_r32_20260922_D/manifest.json`; report regeneration is offline.
- Resume accepts exact fingerprint-matched qualified results. Unclosed intervals, corrupt evidence, and earlier unqualified attempts require review.
