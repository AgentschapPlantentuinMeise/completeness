#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Completeness mapping per 1×1 km EEA cell for plant occurrence data (incidence-based).
Implements Chao et al. (2020) sample completeness estimators for incidence data (q=0 and q=1),
and bootstraps confidence intervals by resampling sampling units within each cell.

Inputs:
- Occurrence file (TSV or CSV): must include eeacellcode, specieskey (or species), yearmonthday
- Flanders boundary polygon (GeoJSON or any vector format readable by GeoPandas)

Outputs:
- GeoPackage with 1 km cell polygons and completeness + CI fields
- Maps:
    1) C1_hat point estimate
    2) C1_lo 95% lower bound
    3) C1_width 95% CI width
"""

import re
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.patches as mpatches
from shapely.geometry import box
import matplotlib.pyplot as plt
from pathlib import Path


# ----------------------------
# Paths (adapt to your setup)
# ----------------------------
BASE = Path("/home/qgroom/Documents/Python Scripts/completeness/")

csv_path = BASE / "onekdata" / "0061226-250827131500795.csv"        # occurrence data
boundary_path = BASE / "borders" / "flanders.geojson"               # Flanders polygon
out_gpkg = BASE / "outputs" / "flanders_completeness_1km_boot.gpkg" # output
out_gpkg.parent.mkdir(parents=True, exist_ok=True)

# Output folder
fig_dir = BASE / "outputs" / "figures"
fig_dir.mkdir(parents=True, exist_ok=True)

out_png = fig_dir / "diagnostic_uncertainty_driven_by_effort.png"

# ---------------------------
# Parameters
# ---------------------------
min_T = 10          # minimum number of sampling units (dates) per cell to keep
n_boot = 100         # bootstrap replications per cell (50 = quick; 100 = paper-like)
seed = 42           # RNG seed
reliable_threshold = 0.80  # choose 0.8 or 0.9 depending on how strict you want to be

# ---------------------------
# Read occurrence data
# ---------------------------
df = pd.read_csv(csv_path, sep="\t", dtype={"eeacellcode": "string"})
df["yearmonthday"] = pd.to_datetime(df["yearmonthday"], errors="coerce")

required = {"eeacellcode", "yearmonthday"}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Choose species identifier
if "specieskey" in df.columns:
    species_col = "specieskey"
elif "species" in df.columns:
    species_col = "species"
else:
    raise ValueError("Need either 'specieskey' or 'species' column in the input data.")

# Drop incomplete records
df = df.dropna(subset=["eeacellcode", species_col, "yearmonthday"]).copy()

# Keep only 1 km cells
df = df[df["eeacellcode"].str.startswith("1km", na=False)].copy()


# ---------------------------
# Define sampling units (incidence replication)
# ---------------------------
# Default: one sampling unit = one day in that cell
df["sample_unit"] = df["yearmonthday"].dt.date

# Deduplicate within sampling unit (presence/absence)
df_inc = df.drop_duplicates(subset=["eeacellcode", species_col, "sample_unit"]).copy()

# ---------------------------
# Save the plots
# ---------------------------

def plot_map(column, title, filename, show=False, dpi=300, vmin=None, vmax=None):
    fig, ax = plt.subplots(figsize=(9, 9))

    flanders_3035.boundary.plot(ax=ax, linewidth=0.6)

    gdf_clip.plot(
        column=column,
        ax=ax,
        legend=True,
        vmin=vmin,
        vmax=vmax
    )

    ax.set_axis_off()
    ax.set_title(title, fontsize=14)

    out_file = fig_dir / filename
    out_file.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(out_file, dpi=dpi, bbox_inches="tight")
    print("Saved:", out_file)

    if show:
        plt.show()
    plt.close(fig)


# ---------------------------
# Helper functions (Chao incidence estimators)
# ---------------------------
def B_hat(Q1, Q2, T):
    """Bias-correction term for incidence coverage completeness (Table 1 note in Chao et al. 2020)."""
    if Q2 > 0:
        return (2 * Q2) / ((T - 1) * Q1 + 2 * Q2)
    if Q1 > 0:
        return 2 / ((T - 1) * (Q1 - 1) + 2)
    return 1.0

def chao2(Sobs, Q1, Q2, T):
    """Classic Chao2 incidence richness estimator."""
    if Q2 > 0:
        return Sobs + ((T - 1) / T) * (Q1 * Q1) / (2 * Q2)
    if Q1 > 1:
        return Sobs + ((T - 1) / T) * (Q1 * (Q1 - 1)) / 2
    return float(Sobs)

def completeness_from_Yi(Yi_vals, T):
    """Compute C1_hat and C0_hat from incidence frequencies Yi for one cell."""
    Yi_vals = np.asarray(Yi_vals, dtype=int)
    Sobs = Yi_vals.size
    Q1 = int(np.sum(Yi_vals == 1))
    Q2 = int(np.sum(Yi_vals == 2))
    U  = int(np.sum(Yi_vals))

    if U == 0 or T <= 0:
        return np.nan, np.nan, Sobs, Q1, Q2, U

    B = B_hat(Q1, Q2, T)
    C1 = 1 - (Q1 / U) * (1 - B)

    S_chao2 = chao2(Sobs, Q1, Q2, T)
    C0 = Sobs / S_chao2 if S_chao2 > 0 else np.nan

    return C1, C0, Sobs, Q1, Q2, U


# ---------------------------
# Point estimates per cell (fast)
# ---------------------------
yi = (df_inc.groupby(["eeacellcode", species_col])["sample_unit"]
            .nunique()
            .rename("Yi")
            .reset_index())

T_per_cell = (df_inc.groupby("eeacellcode")["sample_unit"]
                   .nunique()
                   .rename("T"))

def summarize_cell(g):
    Yi_vals = g["Yi"].to_numpy()
    Sobs = len(Yi_vals)
    Q1 = int(np.sum(Yi_vals == 1))
    Q2 = int(np.sum(Yi_vals == 2))
    U  = int(np.sum(Yi_vals))
    return pd.Series({"Sobs": Sobs, "Q1": Q1, "Q2": Q2, "U": U})

summary = (yi.groupby("eeacellcode")
             .apply(summarize_cell)
             .join(T_per_cell, how="left")
             .reset_index())

summary["Bhat"] = summary.apply(lambda r: B_hat(r["Q1"], r["Q2"], r["T"]), axis=1)
summary["C1_hat"] = 1 - (summary["Q1"] / summary["U"]) * (1 - summary["Bhat"])
summary.loc[summary["U"] == 0, "C1_hat"] = np.nan

summary["S_chao2"] = summary.apply(lambda r: chao2(r["Sobs"], r["Q1"], r["Q2"], r["T"]), axis=1)
summary["C0_hat"] = summary["Sobs"] / summary["S_chao2"]

summary.loc[summary["T"] < min_T, ["C1_hat", "C0_hat"]] = np.nan


# ---------------------------
# Bootstrap confidence intervals per cell
# ---------------------------
def bootstrap_cell(df_cell_inc, n_boot=100, seed=1):
    """
    Nonparametric bootstrap for one cell:
    resample sampling units with replacement, recompute completeness.
    df_cell_inc must have columns: species_col, sample_unit (deduplicated).
    """
    rng = np.random.default_rng(seed)

    units = []
    for su, sub in df_cell_inc.groupby("sample_unit"):
        units.append(sub[species_col].to_numpy())

    T = len(units)
    if T == 0:
        return None

    # observed
    Yi_obs = pd.Series(np.concatenate(units)).value_counts().to_numpy()
    C1_obs, C0_obs, Sobs, Q1, Q2, U = completeness_from_Yi(Yi_obs, T)

    C1_boot = np.empty(n_boot, dtype=float)
    C0_boot = np.empty(n_boot, dtype=float)

    for b in range(n_boot):
        idx = rng.integers(0, T, size=T)
        boot_species = np.concatenate([units[i] for i in idx])
        Yi_boot = pd.Series(boot_species).value_counts().to_numpy()
        C1b, C0b, *_ = completeness_from_Yi(Yi_boot, T)
        C1_boot[b] = C1b
        C0_boot[b] = C0b

    return {
        "T": T,
        "C1_hat_boot": C1_obs,
        "C0_hat_boot": C0_obs,
        "C1_lo": np.nanquantile(C1_boot, 0.025),
        "C1_hi": np.nanquantile(C1_boot, 0.975),
        "C0_lo": np.nanquantile(C0_boot, 0.025),
        "C0_hi": np.nanquantile(C0_boot, 0.975),
    }

print(f"Bootstrapping {n_boot} reps per cell...")
rows = []
for cell, sub in df_inc.groupby("eeacellcode"):
    # Skip low replication cells early (saves time)
    T_cell = sub["sample_unit"].nunique()
    if T_cell < min_T:
        continue

    res = bootstrap_cell(sub, n_boot=n_boot, seed=seed)
    if res is None:
        continue
    res["eeacellcode"] = cell
    rows.append(res)

boot = pd.DataFrame(rows)

# Merge bootstrap results into summary table
summary = summary.merge(
    boot[["eeacellcode", "C1_lo", "C1_hi", "C0_lo", "C0_hi"]],
    on="eeacellcode",
    how="left"
)

# Derived uncertainty metrics
summary["C1_width"] = summary["C1_hi"] - summary["C1_lo"]
summary["C0_width"] = summary["C0_hi"] - summary["C0_lo"]


# ---------------------------
# Generate polygons from eeacellcode
# ---------------------------
EEA_RE = re.compile(r"^(?P<size_km>\d+)\s*kmE(?P<E>\d+)N(?P<N>\d+)$", flags=re.IGNORECASE)

def eea_polygon(eeacellcode: str):
    """Return a shapely polygon in EPSG:3035 from an EEA code like 1kmE3844N3146."""
    if not isinstance(eeacellcode, str):
        return None
    m = EEA_RE.match(eeacellcode.strip())
    if not m:
        return None

    size_km = int(m.group("size_km"))
    e_idx = int(m.group("E"))
    n_idx = int(m.group("N"))

    cell_m = size_km * 1000
    minx = e_idx * 1000
    miny = n_idx * 1000
    maxx = minx + cell_m
    maxy = miny + cell_m
    return box(minx, miny, maxx, maxy)

summary["geometry"] = summary["eeacellcode"].apply(eea_polygon)
gdf_cells = gpd.GeoDataFrame(summary, geometry="geometry", crs="EPSG:3035")
gdf_cells = gdf_cells[~gdf_cells["geometry"].isna()].copy()


# ---------------------------
# Read and clip to Flanders boundary
# ---------------------------
flanders = gpd.read_file(boundary_path)

if flanders.crs is None:
    flanders = flanders.set_crs("EPSG:4326")

flanders_3035 = flanders.to_crs(gdf_cells.crs)
gdf_clip = gpd.clip(gdf_cells, flanders_3035)


# ---------------------------
# Export
# ---------------------------
gdf_clip.to_file(out_gpkg, layer="completeness_1km_bootstrap", driver="GPKG")
print("Saved GeoPackage:", out_gpkg)


# ---------------------------
# Plot maps
# ---------------------------

plot_map("C1_hat",   f"Survey completeness per 1×1 km cell (C1̂; min T={min_T})",
         "C1_hat.png", vmin=0, vmax=1)

plot_map("C1_lo",    f"C1̂ 95% lower bound (min T={min_T})",
         "C1_lower.png", vmin=0, vmax=1)

plot_map("C1_width", f"Uncertainty in C1̂ (95% CI width; min T={min_T})",
         "C1_width.png")

from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Output path
fig_dir = BASE / "outputs" / "figures"
fig_dir.mkdir(parents=True, exist_ok=True)
out_png = fig_dir / "completeness_1km_stacked_4panel.png"

# Make 4 stacked axes
fig, axes = plt.subplots(
    4, 1,
    figsize=(8, 20),
    constrained_layout=True
)

# Shared scale for completeness panels (A & B)
norm_c = Normalize(vmin=0, vmax=1)
cmap_c = "viridis"

gdf_clip["reliable_complete"] = (gdf_clip["C1_lo"] >= reliable_threshold)

# ---------------------------
# A) C1_hat
# ---------------------------
flanders_3035.boundary.plot(ax=axes[0], linewidth=0.6)
gdf_clip.plot(column="C1_hat", ax=axes[0], cmap=cmap_c, norm=norm_c)
axes[0].set_title("A) Coverage completeness (C1̂)")
axes[0].set_axis_off()

# ---------------------------
# B) C1_lo
# ---------------------------
flanders_3035.boundary.plot(ax=axes[1], linewidth=0.6)
gdf_clip.plot(column="C1_lo", ax=axes[1], cmap=cmap_c, norm=norm_c)
axes[1].set_title("B) C1̂ 95% lower bound")
axes[1].set_axis_off()

# ---------------------------
# C) C1_width (uncertainty)
# ---------------------------
flanders_3035.boundary.plot(ax=axes[2], linewidth=0.6)
gdf_clip.plot(column="C1_width", ax=axes[2], cmap="magma", legend=True)
axes[2].set_title("C) Uncertainty (95% CI width)")
axes[2].set_axis_off()

# ---------------------------
# D) Binary reliable completeness mask
# ---------------------------
flanders_3035.boundary.plot(ax=axes[3], linewidth=0.6)

# Plot mask: True vs False
mask_gdf_true = gdf_clip[gdf_clip["reliable_complete"] == True]
mask_gdf_false = gdf_clip[gdf_clip["reliable_complete"] == False]

# Choose colours (simple black/white or green/grey; here neutral)
mask_gdf_false.plot(ax=axes[3], color="lightgrey", linewidth=0)
mask_gdf_true.plot(ax=axes[3], color="black", linewidth=0)

axes[3].set_title(f"D) Reliably complete (C1_lo ≥ {reliable_threshold:.2f})")
axes[3].set_axis_off()

# Add legend for mask
legend_patches = [
    mpatches.Patch(color="black", label="Reliable"),
    mpatches.Patch(color="lightgrey", label="Not reliable")
]
axes[3].legend(handles=legend_patches, loc="lower left", frameon=True)

# ---------------------------
# Shared colour bar for completeness (A & B)
# ---------------------------
sm = ScalarMappable(norm=norm_c, cmap=cmap_c)
sm._A = []
cbar = fig.colorbar(
    sm,
    ax=axes[:2],
    orientation="horizontal",
    fraction=0.04,
    pad=0.02
)
cbar.set_label("Sample coverage (C1̂)")

# Save
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", out_png)



# Table export: no geometry, just codes + metrics
out_table = BASE / "outputs" / "cell_completeness_bootstrap_1km.tsv"

cols = [
    "eeacellcode",
    "T", "Sobs", "Q1", "Q2", "U",
    "C1_hat", "C1_lo", "C1_hi", "C1_width",
    "C0_hat", "C0_lo", "C0_hi", "C0_width",
    "reliable_complete"
]

# Some cells may not have C0 fields if you removed them; use intersection
cols = [c for c in cols if c in gdf_clip.columns]

gdf_clip[cols].to_csv(out_table, sep="\t", index=False)
print("Saved table:", out_table)

pct = 100 * gdf_clip["reliable_complete"].mean()
print(f"{pct:.1f}% of sampled 1 km cells are reliably complete at threshold {reliable_threshold}")

############################################################################

out_png = fig_dir / "diagnostic_uncertainty_driven_by_effort.png"

plot_df = gdf_clip[["C1_hat", "C1_width", "T"]].dropna().copy()
plot_df["logT"] = np.log10(plot_df["T"])

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

# Panel 1: C1_hat vs CI width
sc = axes[0].scatter(plot_df["C1_hat"], plot_df["C1_width"], c=plot_df["logT"], s=10, alpha=0.6)
axes[0].set_xlabel("Coverage completeness (C1̂)")
axes[0].set_ylabel("CI width (C1̂)")
axes[0].set_title("Completeness vs uncertainty")

# Panel 2: T vs CI width
axes[1].scatter(plot_df["T"], plot_df["C1_width"], s=10, alpha=0.4)
axes[1].set_xscale("log")
axes[1].set_xlabel("Sampling effort: T (log scale)")
axes[1].set_ylabel("CI width (C1̂)")
axes[1].set_title("Uncertainty decreases with effort")

# Shared colourbar for panel 1
cbar = fig.colorbar(sc, ax=axes[0], orientation="vertical", fraction=0.05, pad=0.03)
cbar.set_label("log10(T)")

fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.close(fig)

print("Saved:", out_png)