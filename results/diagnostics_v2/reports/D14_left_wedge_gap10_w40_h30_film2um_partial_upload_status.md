# D14 Partial Angle Results Upload Status

Status date: 2026-06-06

The D14 left-wedge angle calculation was stopped locally by request. The latest local checkpoint contains 110 completed angle rows out of 150 planned rows.

## Completed Coverage

- Ez, 8.1014 um: 15/15 angles
- Ez, 9.0016 um: 15/15 angles
- Ez, 10.0090 um: 15/15 angles
- Ez, 10.9850 um: 15/15 angles
- Ez, 12.3450 um: 15/15 angles
- Hz, 8.1014 um: 15/15 angles
- Hz, 9.0016 um: 15/15 angles
- Hz, 10.0090 um: 5/15 angles
- Hz, 10.9850 um: 0/15 angles
- Hz, 12.3450 um: 0/15 angles

## Local Commit

A local commit was created with the D14 scripts, D14 checkpoint tables, D14 report, geometry figures, normal spectra figure, and partial angle-resolved figures:

- Branch: codex/d14-partial-angle-results
- Commit: 09c904d Add D14 partial angle results

The local `git push` failed because this machine has no usable GitHub HTTPS credential and no usable SSH private key. The remote branch was created through the GitHub App, but the full binary figure set is still only in the local commit until credentials are fixed.

## Files In Local Commit

- scripts/D13_tilt20_wrapped_gap10_film2um_angle_screen.py
- scripts/D14_left_wedge_gap10_film2um_angle_screen.py
- results/diagnostics_v2/tables/D14_left_wedge_gap10_w40_h30_film2um_angle_spectra_checkpoint.csv
- results/diagnostics_v2/tables/D14_left_wedge_gap10_w40_h30_film2um_normal_metrics.csv
- results/diagnostics_v2/tables/D14_left_wedge_gap10_w40_h30_film2um_normal_spectra.csv
- results/diagnostics_v2/tables/D14_left_wedge_gap10_w40_h30_film2um_normal_spectra_checkpoint.csv
- results/diagnostics_v2/reports/D14_left_wedge_gap10_w40_h30_film2um_report.md
- D14 geometry, normal spectra, heatmap, A-vs-theta, polar, and coverage PNG figures under results/diagnostics_v2/figures/

## Resume Command After Credential Fix

```powershell
git push -u origin codex/d14-partial-angle-results
```
