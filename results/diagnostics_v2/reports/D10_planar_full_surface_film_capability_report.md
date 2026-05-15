# D10 Planar Full-Surface Film Capability Screen

Overall result level: **CODE_PASS**

## Purpose

This is a material-capability upper-bound test for `measured_lossy_wall_film` on a flat Ti backplane.  It is not the final laser-processed sample geometry and does not rename the film as TiO2.

## Required Answers

1. Was measured n,k used in one-wavelength Meep simulations?
   Yes.

2. Best strict-band mean absorptance in 8.1014-12.398 um:
   best strict mean A = 0.8184 at t=3 um, Ez.

3. Required thickness for best result:
   See best result above.

4. Is the best thickness experimentally justified?
   Classification is reported per row; t>1 um is only a capability upper bound unless experiments confirm it.

5. Thresholds:
   mean_A >= 0.60: yes; mean_A >= 0.80: yes.

6-8. Route decision:
   High-emissivity base-layer candidate found; next add asymmetric microstructure while preserving absorption.

9. Strict vs extended:
   Strict conclusions use only 8.1014-12.398 um.  12.398-12.962 um is an extended observation region with `TI_SUBSTRATE_MATERIAL_EXTRAPOLATION_WARNING`.

## Quality Notes

- Transmission must satisfy abs(T)<=1e-3 for inclusion in metrics.
- Coverage below 90% of planned strict points is marked `INCOMPLETE_NUMERICAL_COVERAGE`.
- 
- Solver source mode must be `single_wavelength_narrowband`; no broadband endpoint normalization is used.

## Metrics Preview

| film_thickness_um | polarization | mean_A_8p1014_12p398_strict | mean_A_8p1014_12p962_extended | enhancement_over_bare_absolute_strict | high_emissivity_route_classification | fabrication_relevance_classification |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | Ez | 0.0702953 | 0.0702953 | 0 | MATERIAL_ROUTE_FAIL_FOR_HIGH_EMISSIVITY | THIN_MODIFIED_LAYER_RANGE |
| 1 | Ez | 0.502405 | 0.502405 | 0.43211 | MATERIAL_HAS_LOSS_BUT_REQUIRES_RESONANT_ARCHITECTURE | POSSIBLE_THICK_MODIFIED_LAYER_RANGE_NEEDS_EXPERIMENT_CONFIRMATION |
| 3 | Ez | 0.818441 | 0.818441 | 0.748146 | HIGH_EMISSIVITY_BASE_LAYER_CANDIDATE | CAPABILITY_UPPER_BOUND_ONLY_NOT_YET_PROCESS_JUSTIFIED |

## Outputs

- Spectra: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/D10_planar_full_surface_film_spectra.csv`
- Metrics: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/D10_planar_full_surface_film_metrics.csv`
- Resolution check: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/D10_planar_full_surface_film_resolution_check.csv`
- Log: `/Users/luckydog/meep_sim/logs/diagnostics_v2/D10_planar_full_surface_film_capability.log`
