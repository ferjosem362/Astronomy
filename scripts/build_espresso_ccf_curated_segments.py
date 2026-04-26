#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

INP = Path("data_processed/espresso/espresso_ccf_header_rv.csv")
OUT = Path("data_processed/espresso/espresso_eso_ccf_curated_rv.csv")
SEGOUT = Path("data_processed/espresso/espresso_eso_ccf_curated_segments_summary.csv")
BADOUT = Path("data_processed/espresso/espresso_eso_ccf_rejected_suspicious_rv.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INP)

# Normalize numeric fields
for c in ["bjd_approx", "mjd_obs", "rv_mps", "rv_err_mps", "berv_mps",
          "ccf_fwhm_mps", "ccf_bis_span_mps", "ccf_contrast"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# Keep only finite core values
df = df[np.isfinite(df["mjd_obs"]) & np.isfinite(df["rv_mps"]) & np.isfinite(df["rv_err_mps"])].copy()

# Main science cut:
# Reject wrong/contaminated CCF family near -18.9 km/s.
# Keep Alpha Cen B-like CCF RVs near -24.8 to -26.7 km/s.
good = df["rv_mps"] < -24000
bad = df[~good].copy()
cur = df[good].copy()

bad.to_csv(BADOUT, index=False)

def classify(row):
    mjd = row["mjd_obs"]
    mode = str(row.get("ins_mode", ""))
    mask = str(row.get("ccf_mask", ""))

    if mode == "SINGLEUHR" or mjd < 60000:
        return "ESPRESSO_UHR", "espresso_2019_uhr_b_only"

    if 60044 <= mjd <= 60085.999:
        return "ESPRESSO_HR_2023_EARLY", "espresso_2023_early_dace_minus_eso_minus39mps"

    if mjd >= 60086:
        return "ESPRESSO_HR_2023PLUS", "espresso_2023plus_nominal"

    return "ESPRESSO_OTHER", "espresso_other"

tmp = cur.apply(classify, axis=1, result_type="expand")
cur["model_inst"] = tmp[0]
cur["espresso_segment"] = tmp[1]

# Output schema compatible with previous RV modeling scripts.
# IMPORTANT:
# rv_mps is already m/s.
# rv_err_mps is already m/s. DO NOT multiply by 1000.
out = pd.DataFrame({
    "time": cur["bjd_approx"],
    "mjd_obs": cur["mjd_obs"],
    "rv": cur["rv_mps"],
    "rv_err": cur["rv_err_mps"],
    "berv": cur.get("berv_mps", np.nan),
    "ccf_fwhm": cur.get("ccf_fwhm_mps", np.nan),
    "ccf_bispan": cur.get("ccf_bis_span_mps", np.nan),
    "ccf_contrast": cur.get("ccf_contrast", np.nan),
    "filename": cur["filename"],
    "ins_mode": cur.get("ins_mode", ""),
    "mask": cur.get("ccf_mask", ""),
    "espresso_segment": cur["espresso_segment"],
    "model_inst": cur["model_inst"],
})

out["rv_centered_by_model_inst"] = out["rv"] - out.groupby("model_inst")["rv"].transform("median")

out = out.sort_values(["time", "filename"]).reset_index(drop=True)
out.to_csv(OUT, index=False)

seg = out.groupby(["model_inst", "espresso_segment"], dropna=False).agg(
    n=("rv", "size"),
    bjd_min=("time", "min"),
    bjd_max=("time", "max"),
    mjd_min=("mjd_obs", "min"),
    mjd_max=("mjd_obs", "max"),
    rv_median=("rv", "median"),
    rv_std=("rv", "std"),
    err_median=("rv_err", "median"),
    fwhm_median=("ccf_fwhm", "median"),
    contrast_median=("ccf_contrast", "median"),
).reset_index()

seg.to_csv(SEGOUT, index=False)

print("Input rows:", len(df))
print("Kept rows:", len(out))
print("Rejected suspicious rows:", len(bad))
print("Saved:", OUT)
print("Saved:", SEGOUT)
print("Saved:", BADOUT)

print("\nCurated segment summary:")
print(seg.to_string(index=False))

print("\nRejected suspicious summary:")
if len(bad):
    bad["night_mjd"] = np.floor(bad["mjd_obs"]).astype("Int64")
    print(bad.groupby("night_mjd").agg(
        n=("rv_mps", "size"),
        mjd_min=("mjd_obs", "min"),
        mjd_max=("mjd_obs", "max"),
        rv_median=("rv_mps", "median"),
        rv_min=("rv_mps", "min"),
        rv_max=("rv_mps", "max"),
        err_median=("rv_err_mps", "median"),
        fwhm_median=("ccf_fwhm_mps", "median"),
        contrast_median=("ccf_contrast", "median"),
    ).reset_index().head(80).to_string(index=False))
else:
    print("None")

print("\nDONE")
