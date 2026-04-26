#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

CCF_PATH = Path("data_processed/espresso/espresso_ccf_header_rv_full.csv")

DACE_CANDIDATES = [
    Path("data_interim/qc/dace_master_full.csv"),
    Path("data_interim/qc/dace_full_master.csv"),
    Path("data_interim/qc/dace_master_nightly.csv"),
    Path("data_interim/qc/espresso_master_full.csv"),
    Path("data_interim/qc/espresso_master.csv"),
]

OUTDIR = Path("data_processed/espresso")
OUTDIR.mkdir(parents=True, exist_ok=True)

def find_dace():
    for p in DACE_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError("No encontré ningún CSV DACE esperado: " + ", ".join(map(str, DACE_CANDIDATES)))

def norm_num(s):
    return pd.to_numeric(s, errors="coerce")

dace_path = find_dace()
print("CCF:", CCF_PATH)
print("DACE:", dace_path)

ccf = pd.read_csv(CCF_PATH)
dace = pd.read_csv(dace_path)

print("\nCCF shape:", ccf.shape)
print("DACE shape:", dace.shape)
print("\nDACE columns:")
print(list(dace.columns))

# Normalizar CCF
ccf["mjd_ccf"] = norm_num(ccf["mjd_obs"])
ccf["bjd_ccf_approx"] = norm_num(ccf["bjd_approx"])
ccf["rv_ccf_mps"] = norm_num(ccf["rv_mps"])
ccf["rv_err_ccf_mps"] = norm_num(ccf["rv_err_mps"])
ccf["ins_mode_ccf"] = ccf.get("hierarch_eso_ins_mode", "").astype(str)
ccf["mask_ccf"] = ccf.get("hierarch_eso_qc_ccf_mask", "").astype(str)

# Normalizar DACE
if "sub_inst" in dace.columns:
    esp = dace[dace["sub_inst"].astype(str).str.contains("ESPRESSO", case=False, na=False)].copy()
else:
    # fallback: si ya es un CSV solo ESPRESSO
    esp = dace.copy()

esp["bjd_dace"] = norm_num(esp["bjd"]) if "bjd" in esp.columns else np.nan
esp["mjd_dace_approx"] = esp["bjd_dace"] - 2400000.5
esp["rv_dace_mps"] = norm_num(esp["rv"]) if "rv" in esp.columns else np.nan
esp["rv_err_dace_mps"] = norm_num(esp["rv_err"]) if "rv_err" in esp.columns else np.nan

print("\nDACE ESPRESSO rows:", len(esp))
if "sub_inst" in esp.columns:
    print("\nDACE ESPRESSO by sub_inst:")
    print(esp["sub_inst"].value_counts(dropna=False).to_string())

print("\nRanges:")
print("CCF  MJD:", ccf["mjd_ccf"].min(), ccf["mjd_ccf"].max())
print("DACE MJD:", esp["mjd_dace_approx"].min(), esp["mjd_dace_approx"].max())

# Match nearest in time. BJD vs MJD_OBS no será exacto por corrección baricéntrica,
# así que usamos tolerancias amplias y diagnosticamos.
ccf2 = ccf.sort_values("bjd_ccf_approx").reset_index(drop=True).copy()
esp2 = esp.sort_values("bjd_dace").reset_index(drop=True).copy()

match = pd.merge_asof(
    esp2,
    ccf2,
    left_on="bjd_dace",
    right_on="bjd_ccf_approx",
    direction="nearest",
    tolerance=0.05,  # 72 min; luego filtramos/diagnosticamos
    suffixes=("_dace", "_ccf"),
)

match["dt_days"] = match["bjd_dace"] - match["bjd_ccf_approx"]
match["dt_min"] = match["dt_days"] * 24 * 60
match["drv_dace_minus_ccf_mps"] = match["rv_dace_mps"] - match["rv_ccf_mps"]
match["derr_quad_mps"] = np.sqrt(match["rv_err_dace_mps"]**2 + match["rv_err_ccf_mps"]**2)

matched = match[match["rv_ccf_mps"].notna()].copy()

# Versiones con tolerancias más estrictas
matched_10min = matched[matched["dt_min"].abs() <= 10].copy()
matched_30min = matched[matched["dt_min"].abs() <= 30].copy()

out_all = OUTDIR / "espresso_dace_vs_eso_ccf_nearest_all.csv"
out_30 = OUTDIR / "espresso_dace_vs_eso_ccf_nearest_30min.csv"
out_10 = OUTDIR / "espresso_dace_vs_eso_ccf_nearest_10min.csv"

match.to_csv(out_all, index=False)
matched_30min.to_csv(out_30, index=False)
matched_10min.to_csv(out_10, index=False)

print("\nSaved:", out_all)
print("Saved:", out_30)
print("Saved:", out_10)

print("\nMatch counts:")
print("DACE ESPRESSO:", len(esp))
print("Matched within 72 min:", len(matched))
print("Matched within 30 min:", len(matched_30min))
print("Matched within 10 min:", len(matched_10min))

def report(label, x):
    print(f"\n== {label} ==")
    if len(x) == 0:
        print("No rows")
        return
    cols = ["dt_min", "drv_dace_minus_ccf_mps", "rv_dace_mps", "rv_ccf_mps", "rv_err_dace_mps", "rv_err_ccf_mps"]
    print(x[cols].agg(["count", "min", "max", "mean", "median", "std"]).to_string())
    if "sub_inst" in x.columns:
        print("\nBy sub_inst:")
        print(x.groupby("sub_inst")["drv_dace_minus_ccf_mps"].agg(["count", "median", "mean", "std", "min", "max"]).to_string())
    if "ins_mode_ccf" in x.columns:
        print("\nBy CCF ins_mode:")
        print(x.groupby("ins_mode_ccf")["drv_dace_minus_ccf_mps"].agg(["count", "median", "mean", "std", "min", "max"]).to_string())

report("Matched 72 min", matched)
report("Matched 30 min", matched_30min)
report("Matched 10 min", matched_10min)

print("\nPreview 30min:")
show_cols = [
    "bjd_dace", "bjd_ccf_approx", "dt_min",
    "rv_dace_mps", "rv_ccf_mps", "drv_dace_minus_ccf_mps",
    "rv_err_dace_mps", "rv_err_ccf_mps",
    "sub_inst", "dace_drs", "dace_mode",
    "ins_mode_ccf", "mask_ccf", "filename",
]
show_cols = [c for c in show_cols if c in matched_30min.columns]
print(matched_30min[show_cols].head(80).to_string(index=False))

print("\nDONE")
