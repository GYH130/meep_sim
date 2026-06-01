# D11 Scaled Deep Slanted Groove Inner-Wall Film Validation

Overall result level: **CODE_PASS**

## Required Answers

1. Scaled geometry built?
   Yes, if geometry checks are PASS: P=50 um, top/bottom width=20/20 um, depth=15 um, tilt=20 deg, substrate=25 um.

2. Is the 500 nm film only on inner sidewalls and bottom?
   Yes by construction: the script uses `build_inner_wall_film_slanted_groove_geometry(..., coating_mode="sidewalls_and_bottom")`, an outer groove film prism, and an inner air prism that opens at the top surface. No top-land film block is added. The geometry PNGs mark the top Ti land as uncoated.

3. 500 nm improvement over scaled bare groove:
   Ez: not available absolute; Hz: not available absolute; 2D proxy: not available absolute.

4. Does the scaled bare groove improve over flat Ti?
   Ez: bare mean not available vs flat not available; Hz: bare mean not available vs flat not available; proxy: bare mean not available vs flat not available.

5. Enhancement source:
   not enough valid metrics to separate geometry and film contributions.

6. Does it reach mean_A_strict >= 0.60?
   No. Best 500 nm strict mean A = 0.255122 (Hz, res=24, PML=4 um).

7. Does it reach mean_A_strict >= 0.80?
   No.

8. If below 0.60:
   This scaled groove + 500 nm inner-wall film structure is not worth treating as the main high-emissivity route for fine optimization, and should not automatically enter angular directionality optimization.

9. If between 0.60 and 0.80:
   Continue with deep-groove size, film thickness, and coverage optimization, and perform stricter resolution validation.

10. If >=0.80:
   Proceed to directionality testing, but because P=50 um supports multiple diffraction orders, use mode decomposition to analyze order-resolved directionality.

11. 2D scope:
   The current result is a 2D equivalent slanted-groove model, not a true 3D laser-processed slanted-hole array.


## Route Decision

SCALED_INNER_WALL_FILM_ROUTE_FAIL_FOR_HIGH_EMISSIVITY or below the continuation threshold. This scaled groove + 500 nm inner-wall film structure should not be treated as the main high-emissivity route and should not automatically enter angle-resolved directionality optimization.

## Quantitative Scope

- Formal metric: `mean_A_8p1014_12p398_strict`.
- Extended observation metric: `mean_A_8p1014_12p962_extended`, flagged with `TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING`.
- `A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2` is only a 2D polarization-average proxy and is not the non-polarized thermal emissivity of a true 3D slanted-pore array.
- Flux monitors cover the full period, so R/T/A are total powers. This is suitable for total absorption screening, not for directional-emission conclusions.

## Geometry Checks

| check_name | status | details |
| --- | --- | --- |
| depth_smaller_than_substrate | PASS | depth=15, substrate=25 |
| outer_groove_inside_unit_cell | PASS | period=50, vertices=[(-10.0, 0.0), (10.0, 0.0), (15.459553513993036, -15.0), (-4.540446486006965, -15.0)] |
| zero_film_degenerates_to_bare_groove | PASS | D11 uses build_slanted_groove_geometry for bare case and the geometry API degenerates film_thickness_um=0 to that same bare geometry. |
| air_core_not_collapsed_250nm | PASS | inner_top=19.5, inner_bottom=19.5, inner_depth=14.75 |
| inner_air_core_inside_unit_cell_250nm | PASS | period=50, vertices=[(-9.75, 0.0), (9.75, 0.0), (15.118560955426485, -14.75), (-4.381439044573516, -14.75)] |
| air_core_not_collapsed_500nm | PASS | inner_top=19, inner_bottom=19, inner_depth=14.5 |
| inner_air_core_inside_unit_cell_500nm | PASS | period=50, vertices=[(-9.5, 0.0), (9.5, 0.0), (14.777568396859934, -14.5), (-4.222431603140066, -14.5)] |
| top_flat_ti_surface_uncoated | PASS | Only an outer groove-cavity film prism and an inner air prism are used; no top-land film blocks are added. |

## Metrics Preview

| case_name | polarization | mean_A_8p1014_12p398_strict | peak_A_strict | enhancement_over_scaled_bare_absolute | enhancement_over_flat_Ti_absolute | route_decision | numerical_status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scaled_inner_wall_film_500nm | Ez | 0.114266 | 0.151332 | nan | nan | SCALED_INNER_WALL_FILM_ROUTE_FAIL_FOR_HIGH_EMISSIVITY | CODE_PASS |
| scaled_inner_wall_film_500nm | Hz | 0.255122 | 0.359113 | nan | nan | SCALED_INNER_WALL_FILM_ROUTE_FAIL_FOR_HIGH_EMISSIVITY | CODE_PASS |
| scaled_inner_wall_film_500nm | unpolarized_2D_proxy | 0.185238 | 0.255222 | nan | nan | SCALED_INNER_WALL_FILM_ROUTE_FAIL_FOR_HIGH_EMISSIVITY | CODE_PASS |

## Numerical Checks

- Resolution check rows: 1
- PML check rows: 1
- Field snapshot status: NOT_RUN
- Field snapshot note: Field maps are only gated from screen-mode metrics.

## Outputs

- Spectra: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D11_scaled_deep_groove_inner_wall_film_spectra.csv`
- Metrics: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D11_scaled_deep_groove_inner_wall_film_metrics.csv`
- Resolution check: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D11_scaled_deep_groove_resolution_check.csv`
- PML check: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/results/diagnostics_v2/tables/D11_scaled_deep_groove_pml_check.csv`
- Log: `/mnt/c/Users/admin/Meep_Simulation/meep_sim/logs/diagnostics_v2/D11_scaled_deep_groove_inner_wall_film_validation.log`
