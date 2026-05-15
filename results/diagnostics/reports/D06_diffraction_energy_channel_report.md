# D06 Diffraction-Order Energy Channel Diagnostic

## Purpose

Use Meep mode decomposition to test whether periodic grooves increase true
absorption or mainly redistribute reflected energy into non-specular
diffraction orders.

## Run Configuration

- Period: 10 um
- Wavelengths: 8 um
- Angles: 0 deg
- Polarizations: Ez
- Resolution: 32 px/um
- PML / air buffer / substrate: 2 / 4 / 4 um
- Decay threshold: 20 dB

## Physical Assumptions And Normalization

- Source propagates from +y toward -y.
- Raw flux monitor convention is retained in the CSV:
  `R = reflection_flux_raw / abs(input_flux_raw)`,
  `T = -transmission_flux_raw / abs(input_flux_raw)`,
  `A = 1 - R - T`.
- Mode decomposition uses `mp.DiffractedPlanewave(g=[m,0,0])`.
  Ez uses s-polarization and Hz uses p-polarization in the x-y incidence plane.
- The script computes both mode-coefficient direction components and selects the
  one whose propagated-order sum best matches the total reflected flux.

## Pass/Fail Criteria

- Finite fluxes and diffraction-order powers.
- `abs(total_reflected_power_from_orders - R_flux_monitor) <= 0.08`.
- For symmetric grooves at normal incidence,
  `abs(sum(+m)-sum(-m)) <= 0.05`.

Overall status: **PASS**

## Allowed Propagating Orders

| wavelength_um | angle_deg | propagating_orders | angles_deg |
|---:|---:|---|---|
| 8 | 0 | -1 0 1 | m=-1:-53.1, m=0:0.0, m=1:53.1 |

## Numerical Findings

Validated by this run:

- Mean absorptance = 0.1346; mean reflected flux = 0.8654. If absorptance remains low while reflected channels sum near R, energy is leaving primarily as reflection rather than being dissipated.
- Slanted nonspecular reflection mean = 0.4278; symmetric mean = 0.4782.
- Mean |positive-negative diffraction order fraction| = 0.0031; max = 0.0062.
- Order-sum residual range:
  -1.021e-08
  to
  -9.138e-09.

Still hypotheses:

- If slanted grooves show larger non-specular reflection but similar A, the
  structure is behaving more like a direction/angle redistribution element than
  a high-emissivity absorber.
- If positive and negative order powers differ in the slanted case, the angular
  asymmetry seen in D05 is likely connected to diffraction-channel imbalance.

Needs higher-fidelity confirmation:

- 3D roughness, oxide layers, rounded groove walls, and experimentally measured
  optical constants are not included.
- Near-grazing diffraction orders can be numerically delicate and should be
  checked with higher resolution/PML before quantitative claims.

## Required Answers

1. Does the current slanted groove transfer specular reflection into
   non-specular orders?
   Slanted nonspecular reflection mean = 0.4278; symmetric mean = 0.4782.

2. Is lack of absorption gain because energy still mainly leaves as reflection?
   Mean absorptance = 0.1346; mean reflected flux = 0.8654. If absorptance remains low while reflected channels sum near R, energy is leaving primarily as reflection rather than being dissipated.

3. Which orders are allowed at 8, 10, 12 um?
   At normal incidence with period 10 um: 8 um allows m=-1,0,+1; 10 um places
   m=+-1 at grazing; 12 um allows only m=0.  The full simulated table is shown
   above.

4. Is directionality from diffraction-order asymmetry rather than absorption
   enhancement?
   Mean |positive-negative diffraction order fraction| = 0.0031; max = 0.0062.

## Output Files

- Orders CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D06_diffraction_orders.csv`
- Energy summary CSV: `/Users/luckydog/meep_sim/results/diagnostics/tables/D06_energy_channel_summary.csv`
- Order bar plot: `/Users/luckydog/meep_sim/results/diagnostics/figures/D06_diffraction_order_barplots.png`
- Energy channel plot: `/Users/luckydog/meep_sim/results/diagnostics/figures/D06_energy_channel_stacked_bars.png`
