#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  1 19:12:55 2026

@author: qgroom
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch occupancy modelling (imperfect detection) using PyMC.
Option A: fit one single-season occupancy model per year per species.

Detectability is modelled as a function of:
  - effort (log survey_days per site-year)
  - completeness proxy: C1_lower per grid cell (incidence-based sample coverage lower bound)

Inputs (from preprocessing):
  - cell_year_survey_days.tsv
  - species_cell_year_detection_days.tsv
  - species_list.tsv
  - C1_lower_5km.tsv  (eeacellcode -> C1_lower)

Outputs:
  - occupancy_trends_all_species_C1lower.tsv
"""

from pathlib import Path
import pandas as pd
import numpy as np

import pymc as pm
import arviz as az


# ---------------------------
# Paths
# ---------------------------
BASE = Path("/home/qgroom/Documents/Python Scripts/occupancy")
in_dir = BASE / "outputs_occupancy"
out_dir = BASE / "outputs_models"
out_dir.mkdir(parents=True, exist_ok=True)

survey_path  = in_dir / "cell_year_survey_days.tsv"
detect_path  = in_dir / "species_cell_year_detection_days.tsv"
species_path = in_dir / "species_list.tsv"
c1_path      = in_dir / "C1_lower_5km.tsv"

# ---------------------------
# Parameters
# ---------------------------
year_start, year_end = 2015, 2025
years = list(range(year_start, year_end + 1))

# Species filtering (important for speed & model stability)
min_cells = 30       # species must occur in at least this many 5km cells (in whole dataset)
min_years = 5        # species must occur in at least this many years
min_records = 100    # species must have at least this many unique species×cell×day records

# Site filtering per year
min_sites_per_year = 50   # minimum cell-years (sites) to fit a year model
min_survey_days = 3       # minimum survey days per cell-year to keep a site (replicates)

# PyMC sampling settings
draws = 500
tune  = 500
chains = 2
cores  = 2
target_accept = 0.9
seed = 42

DATE_SEP = ","


# ---------------------------
# Load preprocessed tables
# ---------------------------
print("Loading inputs...")
survey = pd.read_csv(survey_path, sep="\t")
detect = pd.read_csv(detect_path, sep="\t")
species_list = pd.read_csv(species_path, sep="\t")

survey["year"] = survey["year"].astype(int)
detect["year"] = detect["year"].astype(int)

print("Survey table:", len(survey))
print("Detect table:", len(detect))
print("Species list:", len(species_list))

# Load completeness (C1_lower) lookup
print("Loading completeness lookup:", c1_path)
c1 = pd.read_csv(c1_path, sep="\t")
if "C1_lower" not in c1.columns:
    raise ValueError("C1_lower_5km.tsv must contain a column named 'C1_lower'")

c1 = c1[["eeacellcode", "C1_lower"]].copy()
c1["C1_lower"] = pd.to_numeric(c1["C1_lower"], errors="coerce")

print("Completeness lookup rows:", len(c1))


# ---------------------------
# Helper: build tidy detection table for one species
# ---------------------------
def build_species_detection_table(specieskey, speciesname):
    """
    Build tidy detection table for one species across all years.
    Output columns:
      eeacellcode, year, date, detected, survey_days_y, C1_lower
    """

    d_sp = detect[detect["specieskey"] == specieskey].copy()

    merged = (
        survey
          .merge(c1, on="eeacellcode", how="left")
          .merge(d_sp, on=["eeacellcode", "year"], how="left")
          .rename(columns={"detect_day_list": "detect_days"})
    )

    merged["detect_days"] = merged["detect_days"].fillna("")
    merged["specieskey"] = specieskey
    merged["species"] = speciesname

    rows = []
    for _, r in merged.iterrows():

        # filter out low-effort sites early
        if int(r["survey_days"]) < min_survey_days:
            continue

        survey_days = r["survey_day_list"].split(DATE_SEP)
        detect_days = set(r["detect_days"].split(DATE_SEP)) if r["detect_days"] else set()

        c1_val = r["C1_lower"]
        # if missing, keep as NaN (we'll impute later inside model)
        c1_val = float(c1_val) if pd.notna(c1_val) else np.nan

        for day in survey_days:
            rows.append({
                "specieskey": specieskey,
                "species": speciesname,
                "eeacellcode": r["eeacellcode"],
                "year": int(r["year"]),
                "date": day,
                "detected": 1 if day in detect_days else 0,
                "survey_days_y": int(r["survey_days"]),
                "C1_lower": c1_val
            })

    return pd.DataFrame(rows)


# ---------------------------
# Helper: fit occupancy model for one year
# ---------------------------
def fit_one_year(d_year):
    """
    Fit single-season occupancy model for one year for one species.

    Occupancy ψ: intercept only
    Detection p: intercept + effort + C1_lower

    Returns:
      mean_occ, hdi_low, hdi_high, n_sites
    """
    # Site-level table = cell-year
    sites = d_year[["eeacellcode", "year", "C1_lower"]].drop_duplicates().reset_index(drop=True)
    sites["site_id"] = np.arange(len(sites))
    d_year = d_year.merge(sites, on=["eeacellcode", "year", "C1_lower"], how="left")

    site_id = d_year["site_id"].values.astype(int)
    y = d_year["detected"].astype(int).values
    n_sites = len(sites)

    # Effort covariate (log survey days)
    effort = d_year.groupby("site_id")["survey_days_y"].first().values.astype(float)
    effort = np.log1p(effort)
    effort = (effort - effort.mean()) / (effort.std() + 1e-9)

    # Completeness covariate (C1_lower)
    comp = sites["C1_lower"].values.astype(float)

    # Impute missing completeness to median (conservative)
    med = np.nanmedian(comp)
    comp = np.where(np.isnan(comp), med, comp)

    # standardize
    comp = (comp - comp.mean()) / (comp.std() + 1e-9)

    with pm.Model() as model:

        # Occupancy (ψ): intercept-only for now
        beta0 = pm.Normal("beta0", 0, 1)
        psi = pm.Deterministic("psi", pm.math.sigmoid(beta0))  # scalar
        z = pm.Bernoulli("z", p=psi, shape=n_sites)

        # Detection (p): effort + completeness
        alpha0 = pm.Normal("alpha0", 0, 1)
        alpha_eff = pm.Normal("alpha_eff", 0, 1)
        alpha_comp = pm.Normal("alpha_comp", 0, 1)

        logit_p = alpha0 + alpha_eff * effort + alpha_comp * comp
        p = pm.Deterministic("p", pm.math.sigmoid(logit_p))  # vector length n_sites

        # Observation likelihood
        p_obs = z[site_id] * p[site_id]
        y_obs = pm.Bernoulli("y_obs", p=p_obs, observed=y)

        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            cores=cores,
            target_accept=target_accept,
            random_seed=seed,
            progressbar=False
        )

    # Posterior draws for ψ (scalar)
    psi_post = idata.posterior["psi"].values.reshape(-1)  # (chains*draws,)

    mean_occ = psi_post.mean()
    hdi = az.hdi(psi_post, hdi_prob=0.95)

    return float(mean_occ), float(hdi[0]), float(hdi[1]), n_sites


# ---------------------------
# Main batch loop
# ---------------------------
results = []

# filter species for modelling
species_keep = species_list.query(
    "n_cells >= @min_cells and n_years >= @min_years and n_records >= @min_records"
).copy()
species_keep = species_keep.head(2)

print(f"Species to model: {len(species_keep)}")
print(species_keep.head())

for k, row in enumerate(species_keep.itertuples(index=False), 1):

    specieskey = int(row.specieskey)
    speciesname = row.species

    print(f"\n[{k}/{len(species_keep)}] {speciesname} ({specieskey})")

    # Build detection history table once for this species
    d_sp = build_species_detection_table(specieskey, speciesname)

    if d_sp.empty:
        print("  No data after filtering; skipping.")
        continue

    for yr in years:
        d_year = d_sp[d_sp["year"] == yr].copy()
        if d_year.empty:
            continue

        # number of sites (cell-year)
        n_sites = d_year[["eeacellcode", "year"]].drop_duplicates().shape[0]
        if n_sites < min_sites_per_year:
            continue

        try:
            mean_occ, hdi_low, hdi_high, n_sites = fit_one_year(d_year)

            results.append({
                "specieskey": specieskey,
                "species": speciesname,
                "year": yr,
                "mean_occupancy": mean_occ,
                "hdi_low": hdi_low,
                "hdi_high": hdi_high,
                "n_sites": n_sites,
            })

            print(f"  {yr}: ψ={mean_occ:.3f}  95%HDI=({hdi_low:.3f},{hdi_high:.3f})  sites={n_sites}")

        except Exception as e:
            print(f"  {yr}: FAILED ({e})")
            continue


# ---------------------------
# Save results
# ---------------------------
res = pd.DataFrame(results)
out_path = out_dir / "occupancy_trends_all_species_C1lower.tsv"
res.to_csv(out_path, sep="\t", index=False)

print("\nSaved:", out_path)
print("Done.")
