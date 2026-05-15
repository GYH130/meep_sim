# D03 Numerical Convergence Diagnostic Report

## Run Configuration

- wavelength range: 8-13 um
- nfreq: 3
- baseline: resolution=32, pml=2, air=8, substrate=4, decay=40
- convergence tolerance: 0.01

## Pass / Fail Criteria

- all R/T/A and raw flux columns finite
- each row preserves input_flux_raw, reflection_flux_raw, transmission_flux_raw
- max |delta A| in 8-13 um versus the highest-accuracy value in that single-factor scan <= tolerance

## Required Answers

1. Is current default resolution=32 enough for quantitative conclusions? No / not proven by this run; max difference vs highest resolution = 0.08342976642293508.
2. Does current PML=2 um affect 8-13 um results? No significant effect within tolerance; max difference vs highest PML = 0.0.
3. Recommended unified baseline for later diagnostics: use the highest tested settings for parameters that exceed tolerance (resolution=max, pml=max, air_buffer=max, substrate=max, decay=max).
4. Is the low-emissivity conclusion stable to numerical parameters? Yes; slanted mean A range across scans = [0.104938, 0.168062].

## Verified Numerical Conclusions

- resolution flat_Ti resolution=32: mean A=0.073531, max |delta A|=8.502e-04, PASS.
- resolution flat_Ti resolution=48: mean A=0.074241, max |delta A|=0.000e+00, PASS.
- resolution slanted_Ez resolution=32: mean A=0.104938, max |delta A|=7.809e-03, PASS.
- resolution slanted_Ez resolution=48: mean A=0.110018, max |delta A|=0.000e+00, PASS.
- resolution slanted_Hz resolution=32: mean A=0.143051, max |delta A|=8.343e-02, FAIL.
- resolution slanted_Hz resolution=48: mean A=0.168062, max |delta A|=0.000e+00, PASS.
- pml flat_Ti pml_thickness_um=2: mean A=0.073531, max |delta A|=0.000e+00, PASS.
- pml slanted_Ez pml_thickness_um=2: mean A=0.104938, max |delta A|=0.000e+00, PASS.
- pml slanted_Hz pml_thickness_um=2: mean A=0.143051, max |delta A|=0.000e+00, PASS.
- air_buffer flat_Ti air_buffer_um=8: mean A=0.073531, max |delta A|=0.000e+00, PASS.
- air_buffer slanted_Ez air_buffer_um=8: mean A=0.104938, max |delta A|=0.000e+00, PASS.
- air_buffer slanted_Hz air_buffer_um=8: mean A=0.143051, max |delta A|=0.000e+00, PASS.
- substrate flat_Ti substrate_thickness_um=4: mean A=0.073531, max |delta A|=0.000e+00, PASS.
- substrate slanted_Ez substrate_thickness_um=4: mean A=0.104938, max |delta A|=0.000e+00, PASS.
- substrate slanted_Hz substrate_thickness_um=4: mean A=0.143051, max |delta A|=0.000e+00, PASS.
- decay flat_Ti decay_db=40: mean A=0.073531, max |delta A|=0.000e+00, PASS.
- decay slanted_Ez decay_db=40: mean A=0.104938, max |delta A|=0.000e+00, PASS.
- decay slanted_Hz decay_db=40: mean A=0.143051, max |delta A|=0.000e+00, PASS.

## Hypotheses Still Not Proven

- D03 does not prove which physical mode controls absorption; it only tests numerical sensitivity.
- D03 does not include oxidation, rounded sidewalls, multiscale roughness, or 3D finite structures.

## Needs Experiment Or Higher-Fidelity Model

- Ti optical constants near the Rakić upper bound and laser-induced oxide layers need experimental confirmation.
- If any scan fails tolerance in full runs, later diagnostics should use the stricter numerical settings.

## Output Files

- resolution CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D03_resolution_convergence.csv`
- PML CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D03_pml_convergence.csv`
- buffer/substrate/decay CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D03_buffer_substrate_decay_convergence.csv`
- resolution PNG: `/Users/luckydog/meep_sim/results/diagnostics/figures/D03_resolution_spectra.png`
- summary PNG: `/Users/luckydog/meep_sim/results/diagnostics/figures/D03_convergence_metric_summary.png`
- report: `/Users/luckydog/meep_sim/results/diagnostics/reports/D03_numerical_convergence_report.md`
