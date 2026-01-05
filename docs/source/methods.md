# Methods

This script estimates **sample completeness** for plant occurrence data using
**incidence-based diversity estimators**, applied independently within
1 × 1 km EEA grid cells. The approach follows established methods for
incidence data, particularly those described by Chao and colleagues.

## Incidence-based framework

Records are grouped into **sampling units** defined as unique combinations of
grid cell and date (i.e. *cell-day* units). For each grid cell, an
**incidence matrix** is constructed in which:

- rows represent sampling units,
- columns represent species,
- matrix entries are binary (1 = detected, 0 = not detected).

Within a sampling unit, multiple records of the same species are treated as a
single detection.

This framework allows estimation of species richness and diversity while
accounting for uneven sampling effort.

## Incidence frequencies

For each grid cell, species are summarised by their **incidence frequency**,
defined as the number of sampling units in which a species is detected.

Let:

- `Qₖ` be the number of species observed in exactly `k` sampling units,
- `T` be the total number of sampling units in the grid cell.

In particular:

- `Q₁` (uniques) represents species detected in exactly one sampling unit,
- `Q₂` (duplicates) represents species detected in exactly two sampling units.

These quantities form the basis of the completeness estimators.

## Completeness estimation

Completeness is estimated using Hill numbers for two diversity orders:

### Richness-based completeness (q = 0)

For `q = 0`, completeness reflects how complete the observed species list is
with respect to true species richness. It is sensitive to rare species and
is strongly influenced by the number of uniques (`Q₁`).

Observed richness is compared to an incidence-based richness estimator
(Chao-type), and completeness is calculated as the ratio of observed to
estimated total richness.

### Shannon-based completeness (q = 1)

For `q = 1`, completeness is based on Shannon diversity and reflects how well
*typical* species in the community are represented. This measure is less
dominated by very rare species and often stabilises earlier than `q = 0`.

Shannon diversity is estimated from incidence data using bias-corrected
incidence probabilities, and completeness is defined as the ratio of observed
to estimated diversity.

## Bootstrap confidence intervals

Uncertainty in completeness estimates is quantified using **bootstrap
resampling of sampling units within each grid cell**:

1. Sampling units are resampled with replacement.
2. Completeness is recalculated for each bootstrap replicate.
3. Confidence intervals are derived from the empirical distribution of
   bootstrap estimates (typically 95% bounds).

Cells with few sampling units tend to produce wide confidence intervals,
indicating greater uncertainty in completeness estimates.

## Spatial aggregation and output

Completeness metrics and confidence intervals are attached to the geometry of
each 1 km grid cell. The resulting spatial dataset can be used to:

- map spatial variation in sampling completeness,
- identify poorly sampled areas,
- assess the reliability of downstream biodiversity indicators.

## Assumptions and limitations

- **Sampling unit equivalence**: all sampling units are assumed to be
  comparable in effort, although effort is not explicitly modelled.
- **Detection heterogeneity**: differences in detectability among species or
  visits are not modelled separately from sampling effort.
- **Sparse data**: cells with very few sampling units may yield unstable
  estimates and should be interpreted cautiously or filtered using minimum
  effort thresholds.

## References

Chao, A., et al. (2020). *Quantifying sample completeness and comparing
diversities among assemblages*. Ecological Research, 35, 292–314.
