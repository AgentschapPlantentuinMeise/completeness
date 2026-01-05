# Outputs

For the GeoPackage layer(s), define each field:

- C0_hat (or whatever you name it): completeness estimate for q=0 (species richness)
- C1_hat: completeness estimate for q=1 (Shannon diversity-based)
- C1_lo: lower 95% CI (bootstrap percentile or BCa—whatever you used)
- C1_hi
- C1_width = C1_hi - C1_lo
- n_units: number of sampling units in the cell
- S_obs: observed species richness in cell
- Q1, Q2 (if you compute uniques/duplicates) — define precisely in incidence terms (species occurring in exactly 1 or 2 sampling units)

Also: describe the map outputs (PNG/PDF) and what they show.
