# D04 Field and Absorbed Power Diagnostic

## Purpose
Check whether representative Ti grooves localize electromagnetic energy into lossy Ti regions.

## Physical assumptions
- Length unit: um; Meep frequency: f = 1 / wavelength_um.
- 2D periodic, normal incidence, independent Ez/Hz polarizations.
- Ti: Meep built-in Rakić Drude-Lorentz model; wavelengths above 12.398 um are extrapolated.
- Flux convention: R = reflection_flux_raw / abs(input_flux_raw); T = -transmission_flux_raw / abs(input_flux_raw); A_flux = 1 - R - T.
- Absorbed power density follows Meep's absorbed_power_density example: DFT fields are requested with yee_grid=True and evaluated as 2*pi*f*Im(conj(E).D). Hotspot fractions use positive-clipped density to avoid small negative numerical artifacts dominating ratios.

## Pass/fail criteria
- A_flux and A_volume absolute difference <= 1.
- R/T/A and raw flux fields finite.
- Wall/bottom hotspot ratio >= 0.3 is treated as a mechanism flag, not a solver failure.

## Numerical results
| case_name | polarization | wavelength_um | absorptance_flux | absorptance_volume | volume_flux_abs_difference | pass_or_fail |
| --- | --- | --- | --- | --- | --- | --- |
| flat_Ti | Ez | 10 | 0.0717796 | 0.0720045 | 0.00022483 | PASS |
| slanted_groove | Ez | 10 | 0.111265 | 0.111417 | 0.000152829 | PASS |

## Hotspot summary
| case_name | polarization | mean_wall_fraction | mean_bottom_fraction | max_metal_e2_enhancement | max_wall_e2_enhancement | max_bottom_e2_enhancement |
| --- | --- | --- | --- | --- | --- | --- |
| flat_Ti | Ez | 0 | 0 | 0.00731826 | nan | nan |
| slanted_groove | Ez | 0.675749 | 0.0203692 | 0.153696 | 0.153696 | 0.00147716 |

## Required answers
1. Does the current slanted groove produce clear local field enhancement? Not clearly indicated (max metal |E|^2 / air-median |E|^2 = 0.154).
2. Is the enhanced field located in lossy Ti rather than mainly in air? Partly yes (largest mean wall absorbed-power fraction = 0.676).
3. Does Hz more readily form groove-wall or groove-bottom dissipation than Ez? Not proven by this run.
4. Is low spectral absorption explained by no effective localized loss mode? Not fully supported; localized loss exists but may be insufficient in total power.

## Verified conclusions
- Flux/volume consistency overall status: PASS for the executed set.
- Generated field maps and integral tables retain raw flux and normalization provenance.

## Hypotheses
- If wall/bottom absorbed-power fractions remain small in the full run, the simple 2D groove likely fails to create an efficient localized Ti loss channel in 8-13 um.
- Differences between Ez and Hz remain 2D polarization effects, not direct proof of 3D unpolarized emissivity.

## Needs higher-fidelity confirmation
- 3D finite grooves, oxide layers, roughness, rounded sidewalls, and measured Ti optical constants.
- More converged D04 runs if A_volume differs strongly from A_flux in any case.

## Output files
- `/Users/luckydog/meep_sim/results/diagnostics/tables/D04_absorbed_power_integrals.csv`
- `/Users/luckydog/meep_sim/results/diagnostics/tables/D04_hotspot_metrics.csv`
- Figures generated: 8 PNG files under `/Users/luckydog/meep_sim/results/diagnostics/figures`

## Run configuration
```json
{
  "wavelengths_um": [
    "10"
  ],
  "no_d02_peaks": true,
  "cases": [
    "flat_Ti",
    "slanted_groove"
  ],
  "polarizations": [
    "Ez"
  ],
  "resolution": 32,
  "period_um": 10.0,
  "top_width_um": 4.0,
  "bottom_width_um": 4.0,
  "depth_um": 3.0,
  "pml_thickness_um": 2.0,
  "substrate_thickness_um": 4.0,
  "air_buffer_um": 4.0,
  "field_air_above_um": 1.5,
  "fwidth_fraction": 0.08,
  "decay_db": 20.0,
  "wall_band_um": 0.35,
  "bottom_band_um": 0.35,
  "volume_flux_tol": 1.0,
  "hotspot_ratio_threshold": 0.3
}
```
