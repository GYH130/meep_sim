# D00 Flux Sign And Fresnel Validation Report

## Run Configuration

- wavelength range: 8-13 um
- nfreq: 21
- resolution: 32 pixels/um
- pml_thickness_um: 4
- substrate_thickness_um: 8
- air_buffer_um: 8
- source/monitor convention: source propagates toward -y; raw flux columns are preserved

## Physical Assumptions

- Length unit is um and Meep frequency is f = 1 / wavelength_um.
- Normal incidence, 2D Ez polarization; flat isotropic interfaces are polarization independent at normal incidence.
- The lossless dielectric extends through the lower PML to approximate a semi-infinite transmitted medium.
- Flat Ti is compared with a semi-infinite Fresnel interface; wavelengths above the Ti Rakić upper validity bound are warnings, not pass/fail samples.

## Pass / Fail Criteria

- lossless interface tolerance: 0.03
- flat Ti Fresnel tolerance: 0.05
- flat Ti transmission tolerance: 0.02

## Answers

1. Current transmission formula sign: selected definition is `signed`. Lossless max error for raw T = 1.921e+00; for signed T = 6.760e-04.
2. Flat Ti Meep vs Fresnel: `PASS` over valid Ti wavelengths. max |R-R_Fresnel| = 1.655e-03; max |A-A_Fresnel| = 1.655e-03; max |T_selected| = 7.389e-53.
3. Must src/simulation.py be modified before microstructure diagnostics? No further modification is required based on this run. Public T matches selected definition: True.
4. Old results to recompute: 本次验证选择 signed T。凡是在本次修正前按 raw T 写出的旧结果都应重新计算，包括 flat_ti_spectrum.csv、periodic_groove_spectrum.csv、slanted_groove_spectra.csv 以及依赖这些表的图和报告。

## Numerical Verification

- lossless interface: `PASS`; max |R-R_th| = 1.806e-05, max |T-T_th| = 6.760e-04, max |A_selected| = 6.710e-04.
- flat Ti warning samples outside material validity: 2.

## Verified Conclusions

- The lossless n=1.5 interface selects `signed` transmission normalization.
- The CSV files retain `input_flux_raw`, `reflection_flux_raw`, `transmission_flux_raw`, `T_raw`, `T_signed`, `T_selected`, and `A_selected`.
- The flat Ti comparison is PASS within the Rakić validity range sampled by this run.

## Hypotheses Still Not Proven

- Whether slanted grooves fail because they lack an absorbing cavity mode is not answered by D00; D00 only validates the measurement pipeline.
- Whether the structure mainly redistributes directionality rather than total emissivity requires a separate angular-incidence absorptance diagnostic.
- Whether oxidation or multiscale roughness dominates real femtosecond-laser Ti emissivity remains a modeling hypothesis.

## Needs Experiment Or Higher-Fidelity Model

- Ti optical constants near and beyond 12.398 um should be checked against FTIR or ellipsometry for the processed sample.
- Oxide layers, rounded sidewalls, 3D finite grooves, and polarization averaging require follow-up models.

## Output Files

- lossless CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D00_lossless_interface_validation.csv`
- flat Ti CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D00_flat_ti_fresnel_validation.csv`
- lossless PNG: `/Users/luckydog/meep_sim/results/diagnostics/figures/D00_lossless_RT_validation.png`
- flat Ti PNG: `/Users/luckydog/meep_sim/results/diagnostics/figures/D00_flat_ti_fresnel_comparison.png`
- report: `/Users/luckydog/meep_sim/results/diagnostics/reports/D00_flux_sign_validation_report.md`
