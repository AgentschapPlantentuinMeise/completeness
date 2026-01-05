# Overview

`completeness.py` calculates **sample completeness** for plant occurrence data on a **1 × 1 km EEA grid** (per `eeacellcode`). Its main purpose is to quantify **spatial variation in recording coverage**—i.e., which grid cells are likely well-sampled versus where observed species richness is still strongly limited by under-sampling.

The script uses **incidence-based** methods: within each grid cell, records are grouped into **sampling units** (typically visit-like units such as *cell-day* combinations). Each species is treated as **present/absent within each sampling unit**, rather than using raw record counts. This aligns the workflow with established incidence-based richness and completeness estimators (e.g. Chao-style approaches).

## What it estimates

For each 1 km cell, the script estimates sample completeness for two Hill-number orders:

- **q = 0 (richness-based completeness)**  
  Completeness with respect to *species richness*. This is sensitive to rare species and highlights how many species may still be unobserved.

- **q = 1 (Shannon-based completeness)**  
  Completeness with respect to *typical/common species* (Shannon diversity). This is less dominated by singletons and very rare taxa, and often stabilises earlier than q = 0.

These two measures are complementary:
- **q = 0** answers “how complete is the species list?”
- **q = 1** answers “how complete is the representation of the community, beyond just rare species?”

## Uncertainty and confidence intervals

To quantify uncertainty, the script applies **bootstrap resampling of sampling units within each cell**:

1. Sampling units for a cell are resampled with replacement.
2. Completeness metrics are recalculated for each bootstrap replicate.
3. Confidence intervals (e.g. 95% bounds) are derived from the bootstrap distribution.

This yields point estimates and uncertainty metrics (e.g. lower/upper confidence bounds and/or confidence interval width). Cells with few sampling units typically show wider intervals.

## Inputs

The script expects an occurrence table (TSV/CSV) that includes at least:

- `eeacellcode` — 1 km EEA grid cell identifier
- `specieskey` (or `species`) — taxon identifier
- `yearmonthday` — date field used to define sampling units (commonly cell-day)

A boundary polygon (e.g. Flanders) is used to select the spatial extent and/or clip outputs.

## Outputs

Typical outputs include:

- A **GeoPackage** containing 1 km cell geometries with completeness metrics and confidence intervals.
- Map figures summarising:
  - completeness point estimates (e.g. `C1_hat`)
  - lower confidence bounds (e.g. `C1_lo`)
  - confidence interval width (useful as a diagnostic of estimate stability)

These products are intended both for analysis and for reporting spatial sampling bias patterns.

## Interpretation and intended use

Completeness values are most informative for comparing cells under a consistent:
- sampling unit definition (e.g. “cell-day”),
- taxonomic scope,
- and time window (if filtered).

Low completeness may reflect low survey effort, uneven sampling, or (in some cases) genuinely high underlying richness combined with insufficient visits. Because the method is incidence-based, it is particularly suited to atlas-style occurrence datasets with repeated visits per grid cell.

## Key assumptions and limitations

- **Sampling unit definition matters**: “cell-day” assumes records from the same day represent a single visit-like unit.
- **Detectability is not explicitly modelled**: completeness reflects observed incidence patterns and does not separate effort from detection probability.
- **Sparse cells are uncertain**: cells with few sampling units often have wide confidence intervals and may need minimum-effort thresholds.
