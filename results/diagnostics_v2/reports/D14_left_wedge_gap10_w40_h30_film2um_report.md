# D14 Left-Apex Triangular Wedge Angle Screen

## Geometry

- Period: 50 um.
- Top wedge opening: 40 um.
- Top Ti gap: 10 um.
- Depth: 30 um.
- Left apex: x=-30.9191 um, y=-30 um.
- Apex shift magnitude: 10.9191 um from tan(20 deg).
- Film: 2 um normal-offset shell on sidewalls/apex.
- Representation: clipped/wrapped periodic triangular wedge pieces.

## Geometry Checks

| check_name | status | details |
| --- | --- | --- |
| top_gap10_preserved | PASS | period=50, top_width=40, gap=10 |
| left_apex_wraps_periodically | PASS | apex_x=-30.9191 um, half_period=25 um, shift=10.9191 um |
| film_shell_positive | PASS | inner_top=33.9317 um, inner_depth=25.4488 um |
| wall_normal_film_offset | PASS | left_top=2, right_top=2, apex_left=2, apex_right=2 um |
| depth_below_substrate | PASS | depth=30 um, substrate=45 um |
| wrapped_polygon_pieces_exist | PASS | outer_pieces=2, inner_pieces=2 |

## Outputs

- Normal spectra: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D14_left_wedge_gap10_w40_h30_film2um_normal_spectra.csv`
- Normal metrics: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D14_left_wedge_gap10_w40_h30_film2um_normal_metrics.csv`
- Angle spectra: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D14_left_wedge_gap10_w40_h30_film2um_angle_spectra.csv`
- Angle metrics: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D14_left_wedge_gap10_w40_h30_film2um_angle_metrics.csv`
- Geometry figure: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/figures/D14_left_wedge_gap10_w40_h30_film2um_geometry.png`
- Continuous wedge geometry figure: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/figures/D14_left_wedge_gap10_w40_h30_film2um_geometry_continuous.png`
- Log: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/logs/diagnostics_v2/D14_left_wedge_gap10_w40_h30_film2um.log`
