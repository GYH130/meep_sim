# D13 Wrapped Tilt20 Gap10 Film2um Angle Screen

Overall result level: **NUMERICAL_SCREENING**

## Geometry

- Period: 50 um.
- Top/bottom groove width: 40/40 um.
- Top Ti gap: 10 um.
- Depth: 30 um.
- Tilt: 20 deg.
- Bottom shift: 10.9191 um.
- Film: 2 um on sidewalls and bottom.
- Representation: clipped/wrapped periodic polygon pieces, not enlarged-period geometry.

## Geometry Checks

| check_name | status | details |
| --- | --- | --- |
| raw_bottom_shift_requires_wrap | PASS | bottom_shift=10.9191 um; raw bottom right=30.9191; half_period=25 |
| top_gap_is_10um | PASS | period=50, top_width=40, gap=10 |
| inner_air_core_not_collapsed | PASS | inner_top=36, inner_bottom=36, inner_depth=28 |
| wrapped_outer_pieces_inside_cell | PASS | outer_piece_count=2, pieces=[[(-25.0, -13.73738709727311), (-19.080892972013928, -30.0), (-25.0, -30.0)], [(-20.0, 0.0), (20.0, 0.0), (25.0, -13.73738709727311), (25.0, -30.0), (-9.08089297201393, -30.0)]] |
| wrapped_inner_pieces_inside_cell | PASS | inner_piece_count=2, pieces=[[(-25.0, -19.232341936182355), (-21.808833440546334, -28.0), (-25.0, -28.0)], [(-18.0, 0.0), (18.0, 0.0), (25.0, -19.232341936182355), (25.0, -28.0), (-7.8088334405463335, -28.0)]] |

## Normal-Incidence Emissivity

_No rows._

## Angle-Resolved Summary

_No rows._

Best proxy angle: not available.

## Numerical Scope

- A(lambda, theta) is used as the directional emissivity proxy under reciprocity and local thermal equilibrium.
- This is a 2D wrapped periodic slanted-groove model, not a true 3D hole-array model.
- Angle grid: -70 to 70 deg, step approximately 10 deg.
- Resolution=12, Courant=0.05, retry Courant=0.025, decay_db=20.0, fwidth_fraction=0.08.

## Outputs

- Normal spectra: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D13_tilt20_wrapped_gap10_w40_h30_film2um_normal_spectra.csv`
- Normal metrics: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D13_tilt20_wrapped_gap10_w40_h30_film2um_normal_metrics.csv`
- Angle spectra: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D13_tilt20_wrapped_gap10_w40_h30_film2um_angle_spectra.csv`
- Angle metrics: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D13_tilt20_wrapped_gap10_w40_h30_film2um_angle_metrics.csv`
- Geometry figure: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/figures/D13_tilt20_wrapped_gap10_w40_h30_film2um_geometry.png`
- Log: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/logs/diagnostics_v2/D13_tilt20_wrapped_gap10_w40_h30_film2um.log`
