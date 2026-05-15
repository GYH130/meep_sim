# D01 Solver Regression After Flux Fix

## Run Configuration

- wavelength range: 8-13 um
- nfreq: 5
- resolution: 32 pixels/um
- period / width / depth: 10 / 4 / 3 um
- asymmetric tilt: 20 deg
- pml / substrate / air buffer: 2 / 4 / 8 um

## Pass / Fail Criteria

- negative T tolerance: 0.001
- A range tolerance around [0, 1]: 0.001
- flat Fresnel tolerance: 0.05
- rectangular vs tilt=0 tolerance: 0.05

## Required Answers

1. Did the 8-13 um mean absorptance change significantly after the flux fix? No; max |mean_A_new - mean_A_old_rawT| = 2.620e-04.
2. Should the old 03 slanted-groove conclusion be discarded/recomputed? Yes. 旧的 03 表是在符号修正前生成的，应废弃并重算；本 D01 也用 pre_fix_A_from_raw_T 量化了差异。
3. Rectangular groove vs tilt=0 slanted groove degeneracy check: PASS; max |dR| = 0.000e+00, max |dA| = 0.000e+00.

## Verified Numerical Conclusions

- flat_ti: PASS; <A>_8-13 = 0.073222; max_negative_T = 3.653e-16; A range = [6.182e-02, 8.651e-02]; reference = fresnel_flat_ti.
- rectangular_groove: PASS; <A>_8-13 = 0.106312; max_negative_T = 8.695e-04; A range = [4.195e-02, 1.629e-01]; reference = flat_ti.
- slanted_tilt0: PASS; <A>_8-13 = 0.106312; max_negative_T = 8.695e-04; A range = [4.195e-02, 1.629e-01]; reference = rectangular_groove.
- slanted_tilt20: PASS; <A>_8-13 = 0.105314; max_negative_T = 7.612e-04; A range = [4.596e-02, 1.580e-01]; reference = slanted_tilt0.

## Hypotheses Still Not Proven

- D01 does not prove whether the groove lacks a real absorbing cavity mode; it only verifies solver consistency after D00.
- D01 does not determine polarization dependence; Ez-only spectra remain a modeling limitation.
- D01 does not prove whether directionality changes without total emissivity gain; angular-incidence absorptance is still needed.

## Needs Experiment Or Higher-Fidelity Model

- Ti optical constants near 13 um and processed-surface oxidation still need experiment or improved material models.
- Rounded sidewalls, oxide layers, 3D finite grooves, and rough multiscale structures are not covered by this regression.

## Output Files

- spectra CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D01_solver_regression_spectra.csv`
- metrics CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D01_solver_regression_metrics.csv`
- figure: `/Users/luckydog/meep_sim/results/diagnostics/figures/D01_flat_rect_slanted_comparison.png`
- report: `/Users/luckydog/meep_sim/results/diagnostics/reports/D01_solver_regression_report.md`
