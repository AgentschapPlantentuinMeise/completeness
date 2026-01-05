#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  1 20:42:55 2026

@author: qgroom
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch occupancy modelling with imperfect detection (PyMC) — one model per species.

Occupancy trend:
  ψ_t follows a random walk on the logit scale.

Detection model:
  y_it ~ Binomial(K_it, z_it * p_it)
  logit(p_it) = α0 + α1 * effort_it + α2 * C1_lower_site

Inputs:
  - cell_year_survey_days.tsv
  - species_cell_year_detection_days.tsv
  - species_list.tsv
  - C1_lower_5km.tsv   (eeacellcode, C1_lower)

Outputs:
  - occupancy_timeseries_RW_all_species_C1lower.tsv
    posterior mean + 95% HDI of ψ_t for each species × year
"""

from pathlib import Path
import pandas as pd
import numpy as np

import pymc as pm
import arviz as az
import time

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

out_path = out_dir / "occupancy_timeseries_RW_all_species_C1lower.tsv"


# ---------------------------
# Parameters
# ---------------------------
year_start, year_end = 1975, 2025
years = np.arange(year_start, year_end + 1)
n_years = len(years)

# species filtering (adjust for runtime)
min_cells = 30
min_years = 5
min_records = 100

# site-year filtering
min_survey_days = 3

# PyMC sampling settings (increase for final)
draws = 800
tune  = 800
chains = 4
cores  = 10
target_accept = 0.9
seed = 42


# ---------------------------
# Load inputs
# ---------------------------
print("Loading inputs...")
survey = pd.read_csv(survey_path, sep="\t")
detect = pd.read_csv(detect_path, sep="\t")
species_list = pd.read_csv(species_path, sep="\t")

survey["year"] = survey["year"].astype(int)
detect["year"] = detect["year"].astype(int)

print("Survey rows:", len(survey))
print("Detect rows:", len(detect))
print("Species list:", len(species_list))

# Completeness lookup
c1 = pd.read_csv(c1_path, sep="\t")[["eeacellcode", "C1_lower"]].copy()
c1["C1_lower"] = pd.to_numeric(c1["C1_lower"], errors="coerce")

# Base site-year table (effort + completeness)
site_year_base = survey.merge(c1, on="eeacellcode", how="left")

# Analysis window + minimum survey days
site_year_base = site_year_base[
    (site_year_base["year"] >= year_start) &
    (site_year_base["year"] <= year_end) &
    (site_year_base["survey_days"] >= min_survey_days)
].copy()

print("Filtered site-year rows:", len(site_year_base))

# ---------------------------
# Species filter
# ---------------------------
species_keep = species_list.query(
    "n_cells >= @min_cells and n_years >= @min_years and n_records >= @min_records"
).copy()

species_keep = species_keep.head(50)

print(f"Species to model: {len(species_keep)}")


# ---------------------------
# Fit model for one species
# ---------------------------
def fit_species_rw(specieskey: int, speciesname: str) -> pd.DataFrame:
    """
    Fit random-walk occupancy model for one species (marginalised likelihood: no latent z).
    Returns a dataframe with year-level occupancy posterior summaries.
    """
    import pytensor.tensor as pt  # safe to import here if not already at top

    # species detections aggregated to cell-year
    d_sp = detect.loc[detect["specieskey"] == specieskey, ["eeacellcode", "year", "detection_days"]].copy()

    # merge detection days into site-year table
    d = site_year_base.merge(d_sp, on=["eeacellcode", "year"], how="left")
    d["detection_days"] = d["detection_days"].fillna(0).astype(int)

    # define unique sites for this species (and keep C1_lower at site level)
    sites = d[["eeacellcode", "C1_lower"]].drop_duplicates().reset_index(drop=True)
    sites["site_id"] = np.arange(len(sites), dtype=int)
    d = d.merge(sites[["eeacellcode", "site_id"]], on="eeacellcode", how="left")

    # Year index 0..n_years-1
    year_to_t = {y: i for i, y in enumerate(years)}
    d["t"] = d["year"].map(year_to_t).astype(int)

    # Observations
    y = d["detection_days"].to_numpy(dtype=int)
    K = d["survey_days"].to_numpy(dtype=int)
    site_id = d["site_id"].to_numpy(dtype=int)
    t_idx = d["t"].to_numpy(dtype=int)

    n_obs = len(d)
    n_sites = len(sites)

    # Effort covariate (per observation)
    effort = np.log1p(K).astype(float)
    effort = (effort - effort.mean()) / (effort.std() + 1e-9)

    # Completeness covariate (per site, static)
    comp_site = sites["C1_lower"].to_numpy(dtype=float)
    med = np.nanmedian(comp_site)
    comp_site = np.where(np.isnan(comp_site), med, comp_site)
    comp_site = (comp_site - comp_site.mean()) / (comp_site.std() + 1e-9)

    with pm.Model() as model:

        # ---- occupancy: random walk in logit space ----
        sigma_rw = pm.HalfNormal("sigma_rw", 1.0)

        eta0 = pm.Normal("eta0", 0, 1.5)
        eps = pm.Normal("eps", 0, sigma_rw, shape=n_years - 1)

        eta = pm.Deterministic(
            "eta",
            pm.math.concatenate([[eta0], eta0 + pm.math.cumsum(eps)])
        )

        psi = pm.Deterministic("psi", pm.math.sigmoid(eta))  # length n_years
        psi_t = psi[t_idx]  # per observation

        # ---- detection model ----
        alpha0 = pm.Normal("alpha0", 0, 1.5)
        alpha_eff = pm.Normal("alpha_eff", 0, 1.0)
        alpha_comp = pm.Normal("alpha_comp", 0, 1.0)

        logit_p = alpha0 + alpha_eff * effort + alpha_comp * comp_site[site_id]
        p = pm.Deterministic("p", pm.math.sigmoid(logit_p))  # per observation

        # ---- marginalised likelihood (integrate out z) ----
        # y > 0:  log(psi) + log Binomial(y | K, p)
        # y == 0: log((1-psi) + psi*(1-p)^K)
        log_binom = pm.logp(pm.Binomial.dist(n=K, p=p), y)

        eps_safe = 1e-12
        logp_pos = pt.log(psi_t + eps_safe) + log_binom
        logp_zero = pt.log((1 - psi_t) + psi_t * pt.pow((1 - p), K) + eps_safe)

        logp = pt.switch(pt.gt(y, 0), logp_pos, logp_zero)
        pm.Potential("lik", logp.sum())

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
        "n_sites": n_sites,
        "n_obs": n_obs
    })

    return out

# ---------------------------
# Batch loop
# ---------------------------
all_results = []
failed = []

for i, row in enumerate(species_keep.itertuples(index=False), 1):
    t0 = time.time()
    specieskey = int(row.specieskey)
    speciesname = row.species

    print(f"\n[{i}/{len(species_keep)}] {speciesname} ({specieskey})")

    try:
        res = fit_species_rw(specieskey, speciesname)
        print(f"  took {(time.time()-t0)/60:.1f} minutes")
        all_results.append(res)

        print("  OK. ψ range:",
              f"{res['mean_occupancy'].min():.3f}–{res['mean_occupancy'].max():.3f}")

    except Exception as e:
        print("  FAILED:", e)
        failed.append((specieskey, speciesname, str(e)))
        continue


# ---------------------------
# Save results
# ---------------------------
if all_results:
    out = pd.concat(all_results, ignore_index=True)
    out.to_csv(out_path, sep="\t", index=False)
    print("\nSaved:", out_path)
else:
    print("\nNo results to save.")

# Save failures log
if failed:
    fail_path = out_dir / "failed_species_RW.tsv"
    pd.DataFrame(failed, columns=["specieskey", "species", "error"]).to_csv(
        fail_path, sep="\t", index=False
    )
    print("Saved failures:", fail_path)
