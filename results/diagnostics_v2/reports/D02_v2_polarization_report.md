# D02_v2 Polarization Quantitative Diagnostic

Overall result level: **CODE_PASS**

## Purpose

Re-run the Ez/Hz polarization diagnostic using single-wavelength simulations
and diagnostics_v2 quality filters.  This avoids broadband endpoint source
artifacts and excludes Ti Rakić extrapolation points from quantitative metrics.

## Physical Assumptions

- Length unit is um; Meep frequency is `f = 1 / wavelength_um`.
- Source propagates from +y to -y.
- D00 flux convention is used: `R=refl/abs(input)`, `T=-trans/abs(input)`,
  `A=1-R-T`.
- Ti uses Meep built-in Rakić Drude-Lorentz model; wavelengths above
  12.398 um are marked `MATERIAL_EXTRAPOLATION`.
- `A_unpolarized_2D_proxy=(A_Ez+A_Hz)/2` is not a true 3D sample emissivity.

## Pass / Warning Criteria

- Quantitative metrics only use 8.0-12.25 um, material-valid,
  `TRANSMISSION NUMERICAL_PASS`, and `numerical_quality_flag=VALID` points.
- 12.5 and 13.0 um are observation points only.
- Hz enhancement flag uses `delta_mean_A_Hz_minus_Ez > 0.02`;
  it does not use `abs(delta)`.
- smoke runs never become `PHYSICS_READY`.

本次 run_scope=smoke，只能说明脚本、字段和质量标记可运行；不能作为物理结论通过。

## Required Answers

1. In wavelength-integrated valid-band terms, is Hz truly higher than Ez?
   No case exceeded the Hz enhancement tolerance in this run.

2. Are symmetric and slanted grooves clearly different?
   Valid-band mean-A spread across groove cases/polarizations is nan.

3. What is the current bare-Ti 2D non-polarized proxy?
   未生成 Ez/Hz 配对的非偏振二维代理值。

4. Is the result worth entering D03 convergence validation?
   Not yet; first resolve failed/warning quality rows, insufficient samples, or run only smoke scope.

## Controls

- Flat Ti Ez/Hz mean-A difference: nan; tolerance 0.01.
- Slanted +20/-20 normal-incidence mean-A max difference: nan;
  tolerance 0.01.

## Metrics Table

| case_name | polarization | mean_A_8_12p25_valid | peak_A_valid_band | peak_wavelength_valid_band_um | delta_mean_A_Hz_minus_Ez | delta_peak_A_Hz_minus_Ez | hz_enhancement_flag | valid_sample_count | result_level |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| flat_ti | Ez | nan | 0.0717791 | 10 | nan | nan | False | 1 | WARNING |

## Outputs

- Spectra CSV: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/D02_v2_polarization_spectra.csv`
- Metrics CSV: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/D02_v2_polarization_metrics.csv`
- Flat control: `/Users/luckydog/meep_sim/results/diagnostics_v2/figures/D02_v2_flat_control.png`
- Ez/Hz spectra: `/Users/luckydog/meep_sim/results/diagnostics_v2/figures/D02_v2_Ez_Hz_spectra.png`
- Unpolarized proxy: `/Users/luckydog/meep_sim/results/diagnostics_v2/figures/D02_v2_unpolarized_proxy.png`
- Log: `/Users/luckydog/meep_sim/logs/diagnostics_v2/D02_v2_polarization.log`

## Run Configuration

```json
{
  "cases": [
    "flat_ti"
  ],
  "polarizations": [
    "Ez"
  ],
  "wavelengths_um": [
    "10"
  ],
  "include_observation_points": true,
  "period_um": 10.0,
  "top_width_um": 4.0,
  "bottom_width_um": 4.0,
  "depth_um": 3.0,
  "resolution": 32,
  "pml_thickness_um": 2.0,
  "substrate_thickness_um": 4.0,
  "air_buffer_um": 4.0,
  "decay_db": 0.0,
  "fwidth_fraction": 0.06,
  "hz_enhancement_tol": 0.02,
  "flat_pol_tol": 0.01,
  "mirror_abs_tol": 0.01,
  "slanted_difference_tol": 0.02,
  "min_valid_samples": 1,
  "run_scope": "smoke"
}
```
