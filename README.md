# completeness
 Estimate and map survey completeness of plant occurrence data using incidence-based sample coverage (Chao et al. 2020)

# Spatial survey completeness mapping (EEA 1 km grid)

This repository contains Python code to quantify and map spatial variation in survey completeness and uncertainty in unstructured plant occurrence data. Completeness is estimated using incidence-based sample coverage (C1̂) following Chao _et al._ (2020), with bootstrap confidence intervals computed per 1 × 1 km EEA reference grid cell. The workflow produces GIS-ready outputs (GeoPackage + tabular summaries) and report-ready figures.

The approach is designed to support biodiversity indicator development by making sampling bias and detectability limitations explicit and mappable.

Chao, A., Kubota, Y., Zelený, D., Chiu, C., Li, C., Kusumoto, B., Yasuhara, M., Thorn, S., Wei, C., Costello, M.J., & Colwell, R.K. (2020). **Quantifying sample completeness and comparing diversities among assemblages.** *Ecological Research*, **35**, 292–314. https://doi.org/10.1111/1440-1703.12102

---

## Key concepts

### What is being estimated?
This workflow estimates sample coverage completeness for incidence data (C1̂), which can be interpreted as:

> the probability that a new record from a given grid cell belongs to a species that has already been detected in that cell.

This completeness measure is relatively robust at fine spatial grain compared to richness-based completeness, and is appropriate for heterogeneous, opportunistically collected occurrence data.

### Why bootstrap?
A non-parametric bootstrap (resampling sampling units within each cell) is used to quantify uncertainty and produce:
- 95% confidence intervals for completeness
- uncertainty maps (CI width)
- a binary *“reliably complete”* mask based on the lower CI bound

---

## Inputs

### Occurrence data (tabular)
A TSV/CSV file with at least the following columns:

| Column | Description |
|-------|-------------|
| `eeacellcode` | EEA grid code (e.g. `1kmE3844N3146`) |
| `yearmonthday` | Observation date (e.g. `2017-02-10`) |
| `specieskey` (preferred) or `species` | Species identifier |

The script assumes occurrence data are already assigned to the EEA grid (`eeacellcode`).  
Sampling units for incidence data are defined as unique observation dates per grid cell (can be adapted).

### Boundary polygon (vector)
A polygon layer for the area of interest (e.g. Flanders), readable by GeoPandas:
- GeoJSON, Shapefile, GeoPackage, etc.

This is used to clip the grid-cell polygons to the study area.

---

## Outputs

### 1) GeoPackage (`.gpkg`)
A spatial dataset with one row per EEA grid cell, containing:
- cell polygon geometry (EPSG:3035)
- completeness estimates and uncertainty metrics
- reliability mask for downstream analyses

Example fields:
- `C1_hat`, `C1_lo`, `C1_hi`, `C1_width`
- `T` (number of sampling units)
- `reliable_complete` (binary mask based on `C1_lo ≥ threshold`)

### 2) Tabular summary (`.tsv`)
A non-spatial table listing `eeacellcode` and completeness metrics (for joins or QA).

### 3) Figures (`.png`)
Report-ready maps and diagnostics, including:
- completeness (C1̂)
- lower confidence bound
- uncertainty (CI width)
- binary reliability mask
- scatterplot showing completeness–uncertainty relationship driven by sampling effort

---

## Method overview

1. **Aggregate occurrences by 1 km EEA grid cell** (`eeacellcode`).
2. Define **sampling units** within each cell (default: unique dates).
3. Construct incidence data by deduplicating within cell × species × sampling unit.
4. Compute incidence frequencies per cell:
   - `Q1`: species occurring in exactly one sampling unit (“uniques”)
   - `U`: total incidence counts
   - `T`: number of sampling units
5. Estimate **coverage completeness** (C1̂) following Chao et al. (2020).
6. Bootstrap by **resampling sampling units with replacement** within each cell to obtain 95% CIs.
7. Generate **EEA grid cell polygons directly from `eeacellcode`** (EPSG:3035).
8. Clip to boundary polygon and export results.

---

## Installation

### Recommended (Ubuntu / Debian)
This workflow uses GeoPandas and the geospatial stack.

## Data extraction (GBIF → grid-day datacube)

Input occurrence records were extracted from the GBIF occurrence table and aggregated into a *grid-day* datacube.  
Each record was assigned to a **1 × 1 km EEA reference grid cell** using `GBIF_EEARGCODE(1000, lat, lon, coordinateUncertainty)`.  
A sampling unit was defined as a unique **observation date (YYYY-MM-DD)** within each grid cell.

Records were filtered to:
- the study area (Flanders and Brussels; using administrative `levelXgid` fields)
- `occurrenceStatus = PRESENT`
- Years 1975 to 2025
- `hasGeospatialIssues = FALSE`
- excluding `TAXON_MATCH_FUZZY`
- excluding `basisOfRecord` = `FOSSIL_SPECIMEN` and `LIVING_SPECIMEN`
- requiring `speciesKey`, `year`, `month`, `day`, and coordinates

The resulting datacube provides the input for incidence-based completeness estimation (Chao et al. 2020), where sampling units are defined as observation dates within each grid cell.

```
SELECT
  kingdom,
  kingdomkey,
  phylum,
  phylumkey,
  class,
  classkey,
  "order",
  orderkey,
  family,
  familykey,
  genus,
  genuskey,
  species,
  specieskey,

  -- Date as YYYY-MM-DD (sampling unit used later for incidence completeness)
  
  PRINTF('%04d-%02d-%02d', "year", "month", "day") AS yearmonthday,

  -- 1 km EEA reference grid code, accounting for coordinate uncertainty
  
  GBIF_EEARGCODE(
    1000,
    decimallatitude,
    decimallongitude,
    COALESCE(coordinateuncertaintyinmeters, 1000)
  ) AS eeacellcode,

  -- Counts aggregated per taxon × grid cell × day (windowed per partition)
  
  IF(ISNULL(orderkey), NULL,
    SUM(COUNT(*)) OVER (
      PARTITION BY
        orderkey,
        GBIF_EEARGCODE(1000, decimallatitude, decimallongitude, COALESCE(coordinateuncertaintyinmeters, 1000)),
        PRINTF('%04d-%02d-%02d', "year", "month", "day")
    )
  ) AS ordercount,

  IF(ISNULL(familykey), NULL,
    SUM(COUNT(*)) OVER (
      PARTITION BY
        familykey,
        GBIF_EEARGCODE(1000, decimallatitude, decimallongitude, COALESCE(coordinateuncertaintyinmeters, 1000)),
        PRINTF('%04d-%02d-%02d', "year", "month", "day")
    )
  ) AS familycount,

  IF(ISNULL(genuskey), NULL,
    SUM(COUNT(*)) OVER (
      PARTITION BY
        genuskey,
        GBIF_EEARGCODE(1000, decimallatitude, decimallongitude, COALESCE(coordinateuncertaintyinmeters, 1000)),
        PRINTF('%04d-%02d-%02d', "year", "month", "day")
    )
  ) AS genuscount,

  COUNT(*) AS occurrences,

  MIN(GBIF_TEMPORALUNCERTAINTY(eventdate, eventtime)) AS mintemporaluncertainty,
  MIN(COALESCE(coordinateuncertaintyinmeters, 1000)) AS mincoordinateuncertaintyinmeters

FROM occurrence

WHERE
  -- Spatial filter (Flanders)
  
  (
    occurrence.level0gid IN ('BEL.2_1', 'BEL.1.1.1_1')
    OR occurrence.level1gid IN ('BEL.2_1', 'BEL.1.1.1_1')
    OR occurrence.level2gid IN ('BEL.2_1', 'BEL.1.1.1_1')
    OR occurrence.level3gid IN ('BEL.2_1', 'BEL.1.1.1_1')
  )

  -- Presence only
  
  AND occurrence.occurrencestatus = 'PRESENT'

  -- Plant taxon filter (Tracheophyta / Magnoliopsida etc. via taxon keys)
  
  AND (
    occurrence.taxonkey IN ('196', '220')
    OR occurrence.acceptedtaxonkey IN ('196', '220')
    OR occurrence.kingdomkey IN ('196', '220')
    OR occurrence.phylumkey IN ('196', '220')
    OR occurrence.classkey IN ('196', '220')
    OR occurrence.orderkey IN ('196', '220')
    OR occurrence.familykey IN ('196', '220')
    OR occurrence.genuskey IN ('196', '220')
    OR occurrence.specieskey IN ('196', '220')
  )

  -- Time range
  
  AND (occurrence."year" >= 1975 AND occurrence."year" <= 2025)

  -- Data quality filters
  
  AND occurrence.hasgeospatialissues = FALSE
  AND NOT GBIF_STRINGARRAYCONTAINS(occurrence.issue, 'TAXON_MATCH_FUZZY', TRUE)
  AND NOT occurrence.basisofrecord IN ('FOSSIL_SPECIMEN', 'LIVING_SPECIMEN')

  -- Required fields
  
  AND (
    occurrence.specieskey IS NOT NULL
    AND occurrence."year" IS NOT NULL
    AND occurrence."month" IS NOT NULL
    AND occurrence."day" IS NOT NULL
    AND occurrence.hascoordinate = TRUE
  )

GROUP BY
  occurrence.kingdom,
  occurrence.kingdomkey,
  occurrence.phylum,
  occurrence.phylumkey,
  occurrence.class,
  occurrence.classkey,
  occurrence."order",
  occurrence.orderkey,
  occurrence.family,
  occurrence.familykey,
  occurrence.genus,
  occurrence.genuskey,
  occurrence.species,
  occurrence.specieskey,
  PRINTF('%04d-%02d-%02d', "year", "month", "day"),
  GBIF_EEARGCODE(1000, decimallatitude, decimallongitude, COALESCE(coordinateuncertaintyinmeters, 1000));
