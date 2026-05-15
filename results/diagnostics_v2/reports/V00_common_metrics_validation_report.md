# V00 Common Metrics Validation

Overall status: **NUMERICAL_PASS**

## Purpose

Validate diagnostics_v2 postprocessing and single-wavelength solver utilities
without writing to the original diagnostics directory.

## Acceptance Checks

| check_name | status | value | reference_value | message |
|---|---|---:|---:|---|
| wavelength_integrated_average | NUMERICAL_PASS | 0.4 | 0.333333 | trapezoid wavelength average differs from simple mean |
| signed_transmission_definition | NUMERICAL_PASS | -0.042 | -0.042 | downward incidence T=-transmission_flux_raw/abs(input_flux_raw) |
| low_source_snr_flag | NUMERICAL_PASS | 0.0001 | 0.001 | weak broadband endpoint is excluded from quantitative metrics |
| negative_transmission_failure | NUMERICAL_PASS | -0.042 | 0.005 | T=-0.042 is a FAIL for opaque Ti substrate checks |
| single_wavelength_flat_ti_10um | NUMERICAL_PASS | 0.0717791 | 0.0730661 | single-wavelength flat Ti A matches Fresnel and T is opaque |

## Flat Ti 10 um Single-Wavelength Result

- A_meep: 0.0717791
- A_fresnel: 0.0730661
- abs(A_meep - A_fresnel): 0.001287
- T: -2.68774e-16
- quality_flag: VALID

## Output Files

- Table: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/V00_common_metrics_validation.csv`
- Report: `/Users/luckydog/meep_sim/results/diagnostics_v2/reports/V00_common_metrics_validation_report.md`
- Log: `/Users/luckydog/meep_sim/logs/diagnostics_v2/V00_validate_common_metrics.log`
