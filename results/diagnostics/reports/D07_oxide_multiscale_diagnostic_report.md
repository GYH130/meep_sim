# D07 Oxide And Multiscale Structure Diagnostic

## Purpose

Test whether the weak emissivity enhancement of bare Ti slanted grooves may be
caused by missing oxide or composite surface effects.

## Source-Code Changes

- `src.materials.get_tio2_medium()` was added as a replaceable TiO2 interface.
  The current implementation is `placeholder_demo`: lossless, non-dispersive,
  and explicitly not quantitative.
- `src.geometry.build_oxidized_slanted_groove_geometry()` was added for
  `top_film_only` and `conformal_approx` oxide layouts.  `oxide_thickness_um=0`
  returns the existing bare slanted-groove geometry.

These changes are backward compatible with existing `get_ti_medium()`,
`build_slanted_groove_geometry()`, and `run_periodic_2d_metal_spectrum()`.

## Physical Assumptions

- Lengths are in um and Meep frequency is `f = 1 / wavelength_um`.
- Ti is Meep's Rakić Drude-Lorentz model.
- TiO2 uses placeholder_demo with `n=2.4`.
- The TiO2 placeholder is only a sensitivity model. It has no mid-IR phonon
  absorption and no validated wavelength range.
- A uses the D00/D01 flux convention:
  `R=refl/abs(input)`, `T=-trans/abs(input)`, `A=1-R-T`.

## Pass/Fail Criteria

- All R/T/A and raw flux columns are finite.
- Raw flux columns are preserved.
- Oxide thickness 0 degenerates to the corresponding bare result within
  `1e-06` in 8-13 um mean absorptance.

Overall status: **PASS**

## Verified Numerical Findings

- Best sensitivity enhancement: oxide_slanted_groove_conformal_approx Ez t=0.2 um, mean_A=0.1449, enhancement_over_bare=0.0159.
- Flat/conformal comparison is incomplete in this run.

## Hypotheses, Not Final Quantitative Claims

- A higher placeholder-index oxide response would indicate that surface
  dielectric layers can tune optical coupling, not that real TiO2 gives that
  exact absorptance.
- If conformal groove oxide improves more than flat oxide in the placeholder
  sweep, the real sample should be checked for sidewall/bottom oxide coverage.
- The slant angle is still primarily assessed for directionality by D05/D06;
  this D07 normal-incidence sweep only tests total absorptance sensitivity.

## Needs Experiment Or Better Model

- Reliable mid-IR TiO2/TiOx optical constants or a Drude-Lorentz fit.
- Oxide thickness distribution after femtosecond laser processing.
- Cross-section geometry, roughness, and possible nanoparticle/recast layers.
- 3D morphology; the current model is a 2D periodic approximation.

## Required Answers

1. Does oxide significantly improve 8-13 um mean absorptance?
   Best sensitivity enhancement: oxide_slanted_groove_conformal_approx Ez t=0.2 um, mean_A=0.1449, enhancement_over_bare=0.0159.

2. Is improvement mainly from flat oxide or groove coverage?
   Flat/conformal comparison is incomplete in this run.

3. Does slant mainly affect absorption enhancement or directionality?
   In this script, slant is tested only through normal-incidence total
   absorptance. Directionality remains the domain of D05/D06.

4. What should experiments characterize first?
   Prioritize oxide thickness, oxide composition/phase and mid-IR n,k,
   cross-section groove geometry, then multiscale roughness.

## Output Files

- Spectra CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D07_oxide_thickness_sweep.csv`
- Metrics CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D07_oxide_enhancement_metrics.csv`
- Geometry check: `/Users/luckydog/meep_sim/results/diagnostics/figures/D07_oxide_geometry_check.png`
- Ez spectra: `/Users/luckydog/meep_sim/results/diagnostics/figures/D07_oxide_spectral_comparison_Ez.png`
- Hz spectra: `/Users/luckydog/meep_sim/results/diagnostics/figures/D07_oxide_spectral_comparison_Hz.png`
- Mean A plot: `/Users/luckydog/meep_sim/results/diagnostics/figures/D07_mean_A_vs_oxide_thickness.png`
- Best-case |E|^2: `/Users/luckydog/meep_sim/results/diagnostics/figures/D07_best_case_E2.png`
- Best-case |H|^2: `/Users/luckydog/meep_sim/results/diagnostics/figures/D07_best_case_H2.png`
- Best-case absorbed power: `/Users/luckydog/meep_sim/results/diagnostics/figures/D07_best_case_absorbed_power.png`

Best-case field snapshot status: `{"field_snapshot_status": "completed", "field_snapshot_case": "oxide_slanted_groove_conformal_approx", "field_snapshot_pol": "Ez", "field_snapshot_wavelength_um": 8.0, "field_snapshot_oxide_thickness_um": 0.2}`
