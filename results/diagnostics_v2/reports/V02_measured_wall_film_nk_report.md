# V02 Measured Wall-Film n,k Preparation

Overall result level: **CODE_PASS**

## Inputs

- Workbook: `data/raw/measured_lossy_wall_film_nk.xlsx`
- Sheet: `Sheet2`
- Row handling: first two rows skipped; first three columns renamed to `wavelength_um`, `n`, `k`.

## Validation Summary

- Total numeric data points: 217
- Data wavelength range: 8.1014-24.927 um
- Points in 8.1014-12.962 um: 121
- n range: 0.97997-1.99604
- k range: 0.162859-0.631708
- Negative k present: False
- Allowed for this narrowband quantitative wall-film simulation: yes

## Notes

- The material is named `measured_lossy_wall_film`; no chemical identity is inferred from the Excel data.
- No extrapolation is allowed outside the measured wavelength range.
