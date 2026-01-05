# Input data

`completeness.py` operates on occurrence records that can be grouped into
**incidence-based sampling units** within **1 × 1 km EEA grid cells**.
Only a small number of required fields are needed, but their interpretation
is important for correct results.

## Required columns

The input occurrence table (TSV or CSV) must contain at least the following
columns:

| Column name     | Description |
|-----------------|-------------|
| `eeacellcode`   | Identifier of the 1 km EEA grid cell in which the record occurs. |
| `specieskey` (or `species`) | Taxon identifier. This may be a numeric key (e.g. GBIF `speciesKey`) or a species name string. |
| `yearmonthday`  | Date of the record, used to define sampling units (see below). |

Additional columns may be present and are ignored by the script.

## Definition of a sampling unit

Completeness estimation is **incidence-based**, meaning that the fundamental
data unit is *detection / non-detection per sampling unit*, not record counts.

In this workflow, a **sampling unit** is defined as a unique combination of: (eeacellcode, yearmonthday)


All records from the same grid cell on the same date are treated as belonging
to a single visit-like sampling event. Within a sampling unit, each species
is recorded as either:

- `1` — detected at least once during the sampling unit  
- `0` — not detected during the sampling unit

Multiple records of the same species within a sampling unit do **not** increase
its weight.

## Example input structure

| eeacellcode   | yearmonthday | specieskey |
|---------------|--------------|------------|
| 1kmE123N456   | 2021-05-12   | 12345      |
| 1kmE123N456   | 2021-05-12   | 67890      |
| 1kmE123N456   | 2021-06-03   | 12345      |
| 1kmE124N456   | 2021-05-20   | 22222      |

In this example:

- The first two rows form **one sampling unit** (same cell, same date).
- The third row is a **different sampling unit** (same cell, different date).
- Completeness is calculated independently for each grid cell.

## Temporal scope

All records included in the input file are analysed together unless the input
data have been pre-filtered by date.

If completeness is intended to represent a specific time window (e.g.
post-2000 records only), this filtering should be applied **before** running
the script.

## Spatial extent

A boundary polygon (e.g. Flanders) is used to:

- select relevant grid cells, and/or
- clip spatial outputs for mapping.

Only grid cells intersecting the boundary are included in the final outputs.

## Data quality considerations

- **Consistency of dates** is critical, as dates define sampling units.
- Cells with very few sampling units may produce unstable completeness
  estimates and wide confidence intervals.
- Presence-only datasets without repeated visits per cell are poorly suited
  to incidence-based completeness estimation.

## Recommended preprocessing

Before running `completeness.py`, it is recommended to:

- remove obvious spatial or taxonomic errors,
- standardise date formats,
- ensure that the taxonomic scope is consistent across the dataset,
- optionally exclude cells below a minimum number of sampling units.
