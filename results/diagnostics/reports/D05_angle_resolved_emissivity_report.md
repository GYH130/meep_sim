# D05 Angle-Resolved Absorptance / Emissivity Proxy

## Purpose
Compute true angular absorptance by separate oblique-incidence simulations.

## Physical assumptions
- Length unit: um; Meep frequency: f = 1 / wavelength_um.
- Surface normal is y and periodic direction is x.
- kx = frequency*sin(theta); source phase is exp(i*2*pi*kx*x).
- Each wavelength and angle is simulated separately, not via one broadband fixed-k run.
- emissivity_proxy=A only under reciprocity, local thermal equilibrium, and same direction/polarization Kirchhoff correspondence.
- Ti Rakić validity upper wavelength is 12.398 um; 13 um is extrapolated.

## Pass/fail criteria
- All raw flux and R/T/A values finite.
- flat_Ti and symmetric_groove symmetry residual <= 0.03.
- theta=0 D02 comparison warning threshold <= 0.04.

## Required answers
1. Does the current slanted groove improve total emissivity? Not clearly (flat mean=0.06556, symmetric mean=nan, slanted mean=0.09487).
2. Does it show directional asymmetry? Not clearly (largest band-averaged |A(+theta0)-A(-theta0)|=0.002422).
3. Main directionality: polarization=Ez, theta_max=0 deg, A(+/-theta0) ratio=1.03
4. Current structure looks more like a low/moderate-emissivity diagnostic structure.

## Verified Numerical Conclusions
- Finite flux status: PASS.
- Control symmetry failures: 0.
- D02 theta=0 warnings: 0.

## Band-Averaged Metrics
| metric_scope | case_name | case_label | polarization | wavelength_um | theta_of_max_absorptance | max_absorptance | theta0_deg | absorptance_plus_theta0 | absorptance_minus_theta0 | asymmetry_at_theta0 | ratio_at_theta0 | theta0_absorptance | theta0_delta_vs_D02 | mean_A_theta | mean_directionality_ratio | validation_note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| band_average | flat_Ti | Flat Ti | Ez | nan | 0 | 0.0717791 | 30 | 0.062455 | 0.062455 | 1.44329e-15 | 1 | 0.0717791 | nan | 0.065563 | 1 | band average over simulated wavelengths |
| band_average | slanted_groove | Slanted groove, tilt=20 | Ez | nan | 0 | 0.111224 | 30 | 0.0878968 | 0.0854748 | 0.00242197 | 1.02834 | 0.111224 | nan | 0.0948651 | 1.02834 | band average over simulated wavelengths |

## Symmetry Residuals
| case_name | polarization | wavelength_um | max_symmetry_residual |
| --- | --- | --- | --- |
| flat_Ti | Ez | 10 | 1.44329e-15 |
| slanted_groove | Ez | 10 | 0.00242197 |

## Mirror Residuals
_No rows._

## Hypotheses
- If slanted_groove asymmetry exceeds flat/symmetric controls while mean A stays low, the geometry mainly modulates direction rather than raising total emissivity.
- Stronger high-emissivity behavior may require oxide layers, roughness, multiscale features, or different period/depth.

## Needs Higher-Fidelity Confirmation
- Full 5-degree angle grid, 8-13 um wavelength grid, converged Hz resolution, 3D geometry, oxide/roughness, measured Ti optical constants.

## Output Files
- `/Users/luckydog/meep_sim/results/diagnostics/tables/D05_angle_resolved_absorptance.csv`
- `/Users/luckydog/meep_sim/results/diagnostics/tables/D05_directionality_metrics.csv`
- `/Users/luckydog/meep_sim/results/diagnostics/figures/D05_absorptance_heatmap_Ez.png`
- `/Users/luckydog/meep_sim/results/diagnostics/figures/D05_absorptance_heatmap_Hz.png`
- `/Users/luckydog/meep_sim/results/diagnostics/figures/D05_band_averaged_directionality.png`

## Run Configuration
```json
{
  "wavelengths_um": [
    "10"
  ],
  "full_wavelength_grid": false,
  "angle_min_deg": -30.0,
  "angle_max_deg": 30.0,
  "angle_step_deg": 30.0,
  "cases": [
    "flat_Ti",
    "slanted_groove"
  ],
  "include_mirror": false,
  "polarizations": [
    "Ez"
  ],
  "period_um": 10.0,
  "top_width_um": 4.0,
  "bottom_width_um": 4.0,
  "depth_um": 3.0,
  "resolution": 32,
  "pml_thickness_um": 2.0,
  "substrate_thickness_um": 4.0,
  "air_buffer_um": 4.0,
  "fwidth_fraction": 0.06,
  "decay_db": 20.0,
  "theta0_deg": 30.0,
  "symmetry_tol": 0.03,
  "theta0_d02_tol": 0.04
}
```
