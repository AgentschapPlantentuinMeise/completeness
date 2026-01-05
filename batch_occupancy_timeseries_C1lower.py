#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  1 20:32:16 2026

@author: qgroom
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Time-series occupancy modelling (imperfect detection) using PyMC.
One model per species (2000–2025), with a random-walk occupancy trend.

Inputs (from preprocessing):
  - cell_year_survey_days.tsv
  - species_cell_year_detection_days.tsv
  - species_list.tsv
  - C1_lower_5km.tsv   (eeacellcode, C1_lower)

Outputs:
  - occupancy_timeseries_all_species_C1lower.tsv
    with posterior mean and 95% HDI per year per species
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
year_start, year_end = 2000, 2025
years = np.arange(year_start, year_end + 1)
n_years = len(years)

# Species filtering
min_cells = 30
min_years = 5
min_records = 100

# Site-year filtering
min_survey_days = 3

# PyMC sampling settings
draws = 800
tune  = 800
chains = 2
cores  = 2
target_accept = 0.9
seed = 42

# Output
out_path = out_dir / "occupancy_timeseries_all_species_C1lower.tsv"


# ---------------------------
# Load inputs
# ---------------------------
print("Loading inputs...")
survey = pd.read_csv(survey_path, sep="\t")
detect = pd.read_csv(detect_path, sep="\t")
species_list = pd.read_csv(species_path, sep="\t")

survey["year"] = survey["year"].astype(int)
detect["year"] = detect["year"].astype(int)

print("Survey:", len(survey))
print("Detect:", len(detect))
print("Species list:", len(species_list))

# Completeness lookup
c1 = pd.read_csv(c1_path, sep="\t")[["eeacellcode", "C1_lower"]].copy()
c1["C1_lower"] = pd.to_numeric(c1["C1_lower"], errors="coerce")

# Base site-year table: effort + completeness
site_year_base = survey.merge(c1, on="eeacellcode", how="left")

# Keep only site-years in analysis window and with enough survey days
site_year_base = site_year_base[
    (site_year_base["year"] >= year_start) &
    (site_year_base["year"] <= year_end) &
    (site_year_base["survey_days"] >= min_survey_days)
].copy()

print("Site-year base table:", len(site_year_base))


# ---------------------------
# Species filter
# ---------------------------
species_keep = species_list.query(
    "n_cells >= @min_cells and n_years >= @min_years and n_records >= @min_records"
).copy()
species_keep = species_keep.head(2)


print(f"Species to model: {len(species_keep)}")


# ---------------------------
# Fit model for one species
# ---------------------------
def fit_species_timeseries(specieskey, speciesname):
    """
    Fit time-series occupancy model for one species:
    - occupancy ψ_t (random walk on logit scale)
    - detection p_it depends on effort + completeness (C1_lower)
    Returns dataframe with year-level posterior summary.
    """

    # species detections aggregated to cell-year
    d_sp = detect[detect["specieskey"] == specieskey][
        ["eeacellcode", "year", "detection_days"]
    ].copy()

    # merge into site-year table
    d = site_year_base.merge(d_sp, on=["eeacellcode", "year"], how="left")
    d["detection_days"] = d["detection_days"].fillna(0).astype(int)

    # keep only cell-years where cell was surveyed (already filtered)
    # define site id across full dataset for this species
    sites = d[["eeacellcode"]].drop_duplicates().reset_index(drop=True)
    sites["site_id"] = np.arange(len(sites))
    d = d.merge(sites, on="eeacellcode", how="left")

    # year index 0..n_years-1
    d["t"] = d["year"].map({y:i for i,y in enumerate(years)}).astype(int)

    # Observations
    y = d["detection_days"].values.astype(int)
    K = d["survey_days"].values.astype(int)
    site_id = d["site_id"].values.astype(int)
    t_idx = d["t"].values.astype(int)

    n_sites = len(sites)

    # Effort covariate: log1p(K) standardised
    effort = np.log1p(K).astype(float)
    effort = (effort - effort.mean()) / (effort.std() + 1e-9)

    # Completeness covariate: per site (static)
    comp = d.groupby("site_id")["C1_lower"].first().values.astype(float)
    med = np.nanmedian(comp)
    comp = np.where(np.isnan(comp), med, comp)
    comp = (comp - comp.mean()) / (comp.std() + 1e-9)

    with pm.Model() as model:

        # ---- occupancy random-walk in logit space ----
        sigma_rw = pm.Exponential("sigma_rw", 1.0)

        eta0 = pm.Normal("eta0", 0, 1)
        eps = pm.Normal("eps", 0, sigma_rw, shape=n_years)

        eta = pm.Deterministic("eta", pm.math.cumsum(pm.math.concatenate([[eta0], eps[1:]])))
        psi = pm.Deterministic("psi", pm.math.sigmoid(eta))  # length n_years

        # latent occupancy state per site-year
        z = pm.Bernoulli("z", p=psi[t_idx], shape=len(y))

        # ---- detection model ----
        alpha0 = pm.Normal("alpha0", 0, 1)
        alpha_eff = pm.Normal("alpha_eff", 0, 1)
        alpha_comp = pm.Normal("alpha_comp", 0, 1)

        logit_p = alpha0 + alpha_eff * effort + alpha_comp * comp[site_id]
        p = pm.Deterministic("p", pm.math.sigmoid(logit_p))

        # ---- binomial likelihood: detection-days out of survey-days ----
        p_obs = z * p
        y_obs = pm.Binomial("y_obs", n=K, p=p_obs, observed=y)

        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            cores=cores,
            target_accept=target_accept,
            random_seed=seed,
            progressbar=False
        )

    # Summarise ψ_t posterior
    psi_post = idata.posterior["psi"].values.reshape(-1, n_years)

    mean_psi = psi_post.mean(axis=0)
    hdi = az.hdi(psi_post, hdi_prob=0.95)

    out = pd.DataFrame({
        "specieskey": specieskey,
        "species": speciesname,
        "year": years,
        "mean_occupancy": mean_psi,
        "hdi_low": hdi[:, 0],
        "hdi_high": hdi[:, 1],
        "n_sites": n_sites
    })

    return out


# ---------------------------
# Batch loop
# ---------------------------
all_results = []

for i, row in enumerate(species_keep.itertuples(index=False), 1):
    specieskey = int(row.specieskey)
    speciesname = row.species

    print(f"\n[{i}/{len(species_keep)}] {speciesname} ({specieskey})")

    try:
        res = fit_species_timeseries(specieskey, speciesname)
        all_results.append(res)

        print("  done. mean ψ range:",
              f"{res['mean_occupancy'].min():.3f}–{res['mean_occupancy'].max():.3f}")

    except Exception as e:
        print("  FAILED:", e)
        continue


# ---------------------------
# Save output
# ---------------------------
if all_results:
    out = pd.concat(all_results, ignore_index=True)
    out.to_csv(out_path, sep="\t", index=False)
    print("\nSaved:", out_path)
else:
    print("\nNo results to save.")
