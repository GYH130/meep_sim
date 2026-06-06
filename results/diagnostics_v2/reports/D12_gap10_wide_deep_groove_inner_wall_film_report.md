# D12 Gap10 Wide/Deep Groove Inner-Wall Film Screen

Overall result level: **NUMERICAL_SCREENING**

## Required Answers

1. P=50 um, gap=10 um, top/bottom width=40/40 um, h=30 um, tilt=0 deg geometry built?
   Yes. Geometry checks passed. This straight-slot variant keeps the original D12 P=50 um period. Substrate=45 um.

2. Does the tilted bottom remain fully inside the unit cell?
   Yes. See `outer_groove_inside_unit_cell` and `inner_air_core_inside_unit_cell_*` in the geometry checks. For this nominal geometry, the bottom center offset is about 0 um.

3. Is the film only on inner sidewalls and bottom?
   Yes by construction: the script uses `build_inner_wall_film_slanted_groove_geometry(..., coating_mode="sidewalls_and_bottom")`, an outer groove film prism, and an inner air prism that opens at the top surface. No top-land film block is added. The geometry PNGs mark the top Ti land as uncoated.

4. Gap10 500 nm improvement over bare gap10 groove:
   Ez: not available absolute; Hz: not available absolute; 2D proxy: not available absolute.

5. Does the bare gap10 wide/deep groove improve over flat Ti?
   Ez: bare mean not available vs flat not available; Hz: bare mean not available vs flat not available; proxy: bare mean not available vs flat not available.

6. Which film thickness is best among 0.5, 1, and 2 um?
   Best available film row: Best film strict mean A = 0.739132 (gap10_inner_wall_film_2um, Hz, res=12, PML=4 um).

7. Best result source:
   The best available row identifies whether it came from Ez, Hz, or the 2D proxy: Best film strict mean A = 0.739132 (gap10_inner_wall_film_2um, Hz, res=12, PML=4 um).

8. Comparison with old D11 gap=30 references:
   500 nm proxy delta: 0.0712183; 1 um proxy delta: 0.167962; 2 um proxy delta: 0.316298.

9. Enhancement source:
   not enough valid metrics to separate geometry and film contributions.

10. Does it reach mean_A_strict >= 0.60?
   Yes. Best film strict mean A = 0.739132 (gap10_inner_wall_film_2um, Hz, res=12, PML=4 um).

11. Does it reach mean_A_strict >= 0.80?
   No.

12. If below 0.60:
   If a valid geometry later remains below 0.60, then even with gap=10 um, top width=40 um, bottom width=40 um, and depth=30 um, the inner-wall-only coating route is not suitable as the main high-emissivity route. The next step should shift to full-surface high-absorption film plus microstructure directionality control.

13. If between 0.60 and 0.80:
   Continue with deep-groove size, film thickness, and coverage optimization, and perform stricter resolution validation.

14. If >=0.80:
   Proceed to directionality testing, but because P=50 um supports multiple diffraction orders in this band, use mode decomposition to analyze order-resolved directionality.

15. 2D scope:
   The current result is a 2D equivalent straight-groove model, not a true 3D laser-processed hole array.


## Route Decision

ROUTE_WORTH_REFINEMENT_AND_STRUCTURE_OPTIMIZATION. Continue with deep-groove dimension, film-thickness, and coverage optimization plus stricter resolution/PML checks.

## Quantitative Scope

- Formal metric: `mean_A_8p1014_12p398_strict`.
- Extended observation metric: `mean_A_8p1014_12p962_extended`, flagged with `TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING`.
- `A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2` is only a 2D polarization-average proxy and is not the non-polarized thermal emissivity of a true 3D laser-processed hole array.
- Flux monitors cover the full period, so R/T/A are total powers. This is suitable for total absorption screening, not for directional-emission conclusions.

## Geometry Checks

| check_name | status | details |
| --- | --- | --- |
| depth_smaller_than_substrate | PASS | depth=30, substrate=45 |
| outer_groove_inside_unit_cell | PASS | period=50, vertices=[(-20.0, 0.0), (20.0, 0.0), (20.0, -30.0), (-20.0, -30.0)] |
| tilted_bottom_vertices_inside_unit_cell | PASS | bottom_center_shift=0 um, bottom_left=-20 um, bottom_right=20 um, allowed=[-25, 25] |
| zero_film_degenerates_to_bare_groove | PASS | D12 uses build_slanted_groove_geometry for bare case and the geometry API degenerates film_thickness_um=0 to that same bare geometry. |
| air_core_not_collapsed_500nm | PASS | inner_top=39, inner_bottom=39, inner_depth=29.5 |
| inner_air_core_inside_unit_cell_500nm | PASS | period=50, vertices=[(-19.5, 0.0), (19.5, 0.0), (19.5, -29.5), (-19.5, -29.5)] |
| air_core_not_collapsed_1000nm | PASS | inner_top=38, inner_bottom=38, inner_depth=29 |
| inner_air_core_inside_unit_cell_1000nm | PASS | period=50, vertices=[(-19.0, 0.0), (19.0, 0.0), (19.0, -29.0), (-19.0, -29.0)] |
| air_core_not_collapsed_2000nm | PASS | inner_top=36, inner_bottom=36, inner_depth=28 |
| inner_air_core_inside_unit_cell_2000nm | PASS | period=50, vertices=[(-18.0, 0.0), (18.0, 0.0), (18.0, -28.0), (-18.0, -28.0)] |
| gap_is_10um | PASS | gap=10 um |
| period_original_D12_straight_slot | PASS | period=50 um; straight-slot D12 keeps the physical array pitch at P=50 um. |
| top_width_plus_gap_equals_period | PASS | top_width=40, gap=10, period=50 |
| top_flat_ti_surface_uncoated | PASS | Only an outer groove-cavity film prism and an inner air prism are used; no top-land film blocks are added. |

## Metrics Preview

| case_name | polarization | mean_A_8p1014_12p398_strict | peak_A_strict | enhancement_over_gap10_bare_absolute | enhancement_over_gap30_D11_reference_absolute | enhancement_over_flat_Ti_absolute | route_decision | numerical_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gap10_inner_wall_film_500nm | Ez | 0.138184 | 0.202808 | nan | 0.0239185 | nan | GAP10_WIDE_DEEP_GROOVE_ROUTE_FAIL_FOR_HIGH_EMISSIVITY | NUMERICAL_SCREENING |
| gap10_inner_wall_film_500nm | Hz | 0.374727 | 0.450601 | nan | 0.119606 | nan | PARTIAL_ABSORPTION_ENHANCEMENT_BUT_BELOW_CONTINUE_THRESHOLD | NUMERICAL_SCREENING |
| gap10_inner_wall_film_1um | Ez | 0.390675 | 0.585276 | nan | 0.150973 | nan | PARTIAL_ABSORPTION_ENHANCEMENT_BUT_BELOW_CONTINUE_THRESHOLD | NUMERICAL_SCREENING |
| gap10_inner_wall_film_1um | Hz | 0.550093 | 0.655802 | nan | 0.18495 | nan | PARTIAL_ABSORPTION_ENHANCEMENT_BUT_BELOW_CONTINUE_THRESHOLD | NUMERICAL_SCREENING |
| gap10_inner_wall_film_2um | Ez | 0.726493 | 0.794044 | nan | 0.31656 | nan | ROUTE_WORTH_REFINEMENT_AND_STRUCTURE_OPTIMIZATION | NUMERICAL_SCREENING |
| gap10_inner_wall_film_2um | Hz | 0.739132 | 0.808292 | nan | 0.316036 | nan | ROUTE_WORTH_REFINEMENT_AND_STRUCTURE_OPTIMIZATION | NUMERICAL_SCREENING |
| gap10_inner_wall_film_500nm | unpolarized_2D_proxy | 0.256456 | 0.326704 | nan | 0.0712183 | nan | GAP10_WIDE_DEEP_GROOVE_ROUTE_FAIL_FOR_HIGH_EMISSIVITY | NUMERICAL_SCREENING |
| gap10_inner_wall_film_1um | unpolarized_2D_proxy | 0.470384 | 0.620539 | nan | 0.167962 | nan | PARTIAL_ABSORPTION_ENHANCEMENT_BUT_BELOW_CONTINUE_THRESHOLD | NUMERICAL_SCREENING |
| gap10_inner_wall_film_2um | unpolarized_2D_proxy | 0.732813 | 0.801168 | nan | 0.316298 | nan | ROUTE_WORTH_REFINEMENT_AND_STRUCTURE_OPTIMIZATION | NUMERICAL_SCREENING |

## Numerical Checks

- Resolution check rows: 1
- PML check rows: 1
- Field snapshot status: NOT_RUN
- Field snapshot note: Field snapshot skipped after completed spectra/metrics because the D12 straight-slot field snapshot was numerically unstable.

## Outputs

- Spectra: `C:\Users\admin\Meep_Simulation\meep_sim\results\diagnostics_v2\tables\D12_gap10_wide_deep_groove_inner_wall_film_spectra.csv`
- Metrics: `C:\Users\admin\Meep_Simulation\meep_sim\results\diagnostics_v2\tables\D12_gap10_wide_deep_groove_inner_wall_film_metrics.csv`
- Resolution check: `C:\Users\admin\Meep_Simulation\meep_sim\results\diagnostics_v2\tables\D12_gap10_wide_deep_groove_resolution_check.csv`
- PML check: `C:\Users\admin\Meep_Simulation\meep_sim\results\diagnostics_v2\tables\D12_gap10_wide_deep_groove_pml_check.csv`
- Log: `C:\Users\admin\Meep_Simulation\meep_sim\logs\diagnostics_v2\D12_gap10_wide_deep_groove_inner_wall_film.log`
