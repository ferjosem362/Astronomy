#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

CCF_PATH = Path("data_processed/espresso/espresso_ccf_header_rv_full.csv")
DACE_CANDIDATES = [
    Path("data_interim/qc/dace_master_full.csv"),
    Path("data_interim/qc/dace_full_master.csv"),
    Path("data_interim/qc/espresso_master_full.csv"),
    Path("data_interim/qc/espresso_master.csv"),
]

OUTDIR = Path("data_processed/espresso")
OUTDIR.mkdir(parents=True, exist_ok=True)

def find_existing(paths):
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError("No input DACE/ESPRESSO file found")

def num(s):
    return pd.to_numeric(s, errors="coerce")

dace_path = find_existing(DACE_CANDIDATES)

ccf = pd.read_csv(CCF_PATH)
dace = pd.read_csv(dace_path)

# ---------- ESO CCF ----------
ccf["mjd_ccf"] = num(ccf["mjd_obs"])
ccf["rv_ccf_mps"] = num(ccf["rv_mps"])
ccf["rv_err_ccf_mps"] = num(ccf["rv_err_mps"])
ccf["night_mjd"] = np.floor(ccf["mjd_ccf"]).astype("Int64")
ccf["ins_mode_ccf"] = ccf.get("hierarch_eso_ins_mode", "").astype(str)
ccf["mask_ccf"] = ccf.get("hierarch_eso_qc_ccf_mask", "").astype(str)

ccf_hr = ccf[
    ccf["ins_mode_ccf"].eq("SINGLEHR")
    & ccf["mask_ccf"].eq("K2")
    & ccf["rv_ccf_mps"].notna()
    & ccf["mjd_ccf"].notna()
].copy()

ccf_hr = ccf_hr.sort_values(["night_mjd", "mjd_ccf", "filename"]).copy()
ccf_hr["idx_in_night"] = ccf_hr.groupby("night_mjd").cumcount()

# ---------- DACE ----------
if "sub_inst" in dace.columns:
    esp = dace[
        dace["sub_inst"].astype(str).str.contains("ESPRESSO19_SINGLEHR11", case=False, na=False)
    ].copy()
else:
    esp = dace.copy()

esp["bjd_dace"] = num(esp["bjd"])
esp["mjd_dace_approx"] = esp["bjd_dace"] - 2400000.5
esp["rv_dace_mps"] = num(esp["rv"])
esp["rv_err_dace_mps"] = num(esp["rv_err"])
esp["night_mjd"] = np.floor(esp["mjd_dace_approx"]).astype("Int64")

esp = esp[
    esp["bjd_dace"].notna()
    & esp["rv_dace_mps"].notna()
    & esp["night_mjd"].notna()
].copy()

esp = esp.sort_values(["night_mjd", "bjd_dace"]).copy()
esp["idx_in_night"] = esp.groupby("night_mjd").cumcount()

# ---------- group counts ----------
ccf_counts = ccf_hr.groupby("night_mjd").size().rename("n_ccf_night")
esp_counts = esp.groupby("night_mjd").size().rename("n_dace_night")
night_cmp = pd.concat([ccf_counts, esp_counts], axis=1).fillna(0).astype(int)
night_cmp["delta_ccf_minus_dace"] = night_cmp["n_ccf_night"] - night_cmp["n_dace_night"]

# Añadimos counts solo una vez: después del merge, para evitar sufijos raros.
matched = esp.merge(
    ccf_hr,
    on=["night_mjd", "idx_in_night"],
    how="inner",
    suffixes=("_dace", "_ccf"),
)

matched = matched.merge(
    night_cmp.reset_index(),
    on="night_mjd",
    how="left",
)

matched["drv_dace_minus_ccf_mps"] = matched["rv_dace_mps"] - matched["rv_ccf_mps"]
matched["dt_naive_min"] = (
    matched["bjd_dace"] - (matched["mjd_ccf"] + 2400000.5)
) * 1440.0

strict = matched[
    matched["n_ccf_night"].eq(matched["n_dace_night"])
].copy()

out_all = OUTDIR / "espresso_dace_vs_eso_ccf_by_night_index_all.csv"
out_strict = OUTDIR / "espresso_dace_vs_eso_ccf_by_night_index_strict_equal_counts.csv"
matched.to_csv(out_all, index=False)
strict.to_csv(out_strict, index=False)

print("CCF:", CCF_PATH)
print("DACE:", dace_path)
print("CCF HR K2 rows:", len(ccf_hr))
print("DACE ESPRESSO19_SINGLEHR11 rows:", len(esp))
print("Matched by night+index:", len(matched))
print("Strict equal-count matched:", len(strict))
print("Saved:", out_all)
print("Saved:", out_strict)

print("\nNight count comparison:")
print(
    night_cmp.value_counts(
        ["n_ccf_night", "n_dace_night", "delta_ccf_minus_dace"]
    ).sort_index().to_string()
)

print("\nProblem nights preview:")
problem = night_cmp[night_cmp["delta_ccf_minus_dace"] != 0]
if len(problem):
    print(problem.head(100).to_string())
else:
    print("No problem nights")

def report(label, x):
    print(f"\n== {label} ==")
    if len(x) == 0:
        print("No rows")
        return
    base_cols = [
        "drv_dace_minus_ccf_mps",
        "dt_naive_min",
        "rv_dace_mps",
        "rv_ccf_mps",
        "rv_err_dace_mps",
        "rv_err_ccf_mps",
    ]
    base_cols = [c for c in base_cols if c in x.columns]
    print(x[base_cols].agg(["count", "min", "max", "mean", "median", "std"]).to_string())

    print("\nBy n_ccf/n_dace:")
    print(
        x.groupby(["n_ccf_night", "n_dace_night"])["drv_dace_minus_ccf_mps"]
        .agg(["count", "median", "mean", "std", "min", "max"])
        .to_string()
    )

report("ALL night+index", matched)
report("STRICT equal-count nights", strict)

print("\nFirst strict rows:")
cols = [
    "night_mjd", "idx_in_night", "n_ccf_night", "n_dace_night",
    "bjd_dace", "mjd_ccf", "dt_naive_min",
    "rv_dace_mps", "rv_ccf_mps", "drv_dace_minus_ccf_mps",
    "rv_err_dace_mps", "rv_err_ccf_mps",
    "sub_inst", "ins_mode_ccf", "mask_ccf", "filename",
]
cols = [c for c in cols if c in strict.columns]
print(strict[cols].head(80).to_string(index=False))

print("\nDONE")
