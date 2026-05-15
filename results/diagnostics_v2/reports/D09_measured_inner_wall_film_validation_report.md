# D09 Measured Inner-Wall Film Validation

Overall result level: **WARNING**

## Scope

- 2D periodic slanted Ti groove model only; not a final 3D slanted-pore prediction.
- The wall film is named `measured_lossy_wall_film`; no chemical identity is inferred.
- Ti substrate uses the existing Rakić model.
- Film n,k are used only through narrowband single-wavelength conductivity media.
- Flux-derived total absorptance A is the quantitative criterion.

## Required Answers

1. Was the measured n,k table read, interpolated, and used in Meep?
   Yes, for coated cases.

2. What is the actual quantitative evaluation interval?
   8.1014-12.962 um, constrained by the measured film data overlap.

3. Which 200/250/300 nm film improves most?
   Best coated case: t=0.25 um, Ez, mean_A=0.1108, absolute enhancement=0.00664, relative enhancement=0.06375.

4. Absolute and relative enhancement over bare slanted groove:
   Best coated case: t=0.25 um, Ez, mean_A=0.1108, absolute enhancement=0.00664, relative enhancement=0.06375.

5. Is there a mean_A >= 0.50 high-emissivity candidate?
   no.

6. If enhancement remains limited, likely causes include insufficient film n,k loss, limited coated area, weak geometric coupling, 2D-vs-3D mismatch, or missing multiscale roughness.

7. Next steps should be selected from angle-resolved validation, period/depth/width scan, film-thickness optimization, and experimental cross-section/composition validation.

8. The film is not automatically TiO2 unless the Excel chemistry is independently confirmed.

## Quality Notes

- Field snapshot status: NOT_RUN
- Minimum resolution in this run: 32; minimum decay_db: 0.0.
- A screen-mode run is only labeled NUMERICAL_SCREENING when resolution>=64 and decay_db>=60.
- Conductivity film field maps are qualitative only; no D04 Lorentz/Drude absorbed-power formula is used as formal film absorption evidence.
- Points with opaque-substrate transmission FAIL are excluded from band metrics.
- Wavelengths above Ti Rakić validity are explicitly flagged in the spectra table.

## Metrics Preview

| case_name | polarization | film_thickness_um | mean_A_data_covered_band | enhancement_over_bare_absolute | emissivity_candidate_class | enhancement_class |
| --- | --- | --- | --- | --- | --- | --- |
| flat_Ti | Ez | 0 | 0.0702947 | -0.0338555 | LOW_ABSORPTION | LIMITED_ENHANCEMENT |
| bare_slanted_groove | Ez | 0 | 0.10415 | 0 | LOW_ABSORPTION | LIMITED_ENHANCEMENT |
| inner_wall_film_slanted_groove | Ez | 0.25 | 0.11079 | 0.00663992 | LOW_ABSORPTION | LIMITED_ENHANCEMENT |
| flat_Ti | Hz | 0 | 0.0702947 | nan | LOW_ABSORPTION | NOT_QUANTITATIVE |
| bare_slanted_groove | Hz | 0 | nan | nan | NOT_QUANTITATIVE | NOT_QUANTITATIVE |
| inner_wall_film_slanted_groove | Hz | 0.25 | nan | nan | NOT_QUANTITATIVE | NOT_QUANTITATIVE |
| flat_Ti | unpolarized_2D_proxy | 0 | 0.0702947 | nan | LOW_ABSORPTION | nan |
| bare_slanted_groove | unpolarized_2D_proxy | 0 | nan | nan | nan | nan |
| inner_wall_film_slanted_groove | unpolarized_2D_proxy | 0.25 | nan | nan | nan | nan |

## Output Files

- Spectra: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/D09_measured_inner_wall_film_spectra.csv`
- Metrics: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/D09_measured_inner_wall_film_metrics.csv`
- Resolution check: `/Users/luckydog/meep_sim/results/diagnostics_v2/tables/D09_resolution_check.csv`
- Log: `/Users/luckydog/meep_sim/logs/diagnostics_v2/D09_measured_inner_wall_film_validation.log`
