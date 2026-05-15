# D02 Polarization Diagnostic Report

## Run Configuration

- wavelength range: 5-15 um
- nfreq: 5
- resolution: 32 pixels/um
- geometry: P=10, top=4, bottom=4, depth=3 um
- include tilt=30: False

## Pass / Fail Criteria

- all R/T/A and raw flux columns finite
- raw flux columns retained for every wavelength sample
- flat Ti polarization consistency: |mean_A_Hz - mean_A_Ez| <= 0.02
- microstructure Hz enhancement flag: mean_A_Hz - mean_A_Ez > 0.02

## Required Answers

1. Flat Ti Ez/Hz consistency: PASS; delta mean A = -2.028e-07.
2. Does Hz show stronger microstructure absorption than Ez? No clear Hz enhancement by the configured threshold.
3. Does the previous Ez-only 'structure enhancement is not obvious' conclusion still hold? Yes under this 2D diagnostic.

## Verified Numerical Conclusions

- Flat Ti: mean A Ez=0.071781, Hz=0.071781, delta=-2.028e-07; peak A Ez=0.071781, Hz=0.071781.
- Slanted groove, tilt=20: mean A Ez=0.111853, Hz=0.126395, delta=1.454e-02; peak A Ez=0.111853, Hz=0.126395.
- Symmetric groove, tilt=0: mean A Ez=0.114011, Hz=0.129244, delta=1.523e-02; peak A Ez=0.114011, Hz=0.129244.

## Hypotheses Still Not Proven

- A higher Hz response in 2D would not automatically imply higher real non-polarized emission.
- This diagnostic does not include 3D finite grooves, roughness, oxide layers, or angular incidence.
- The result does not identify the microscopic absorption mode; field maps or mode diagnostics are needed.

## Needs Experiment Or Higher-Fidelity Model

- Non-polarized real samples require TE/TM averaging and 3D geometry validation.
- Ti optical constants beyond the Rakić range and laser-induced oxidation need experimental confirmation.

## Output Files

- spectra CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D02_polarization_spectra.csv`
- metrics CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D02_polarization_metrics.csv`
- flat control PNG: `/Users/luckydog/meep_sim/results/diagnostics/figures/D02_flat_Ti_Ez_Hz_control.png`
- structured spectra PNG: `/Users/luckydog/meep_sim/results/diagnostics/figures/D02_slanted_Ez_Hz_spectra.png`
- mean-A bar PNG: `/Users/luckydog/meep_sim/results/diagnostics/figures/D02_mean_A_polarization_comparison.png`
- report: `/Users/luckydog/meep_sim/results/diagnostics/reports/D02_polarization_diagnostic_report.md`
