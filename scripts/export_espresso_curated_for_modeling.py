#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

INP = Path("data_processed/espresso/espresso_eso_ccf_curated_rv.csv")
OUTDIR = Path("data_interim/qc/periodogram_inputs")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INP)

for c in ["time", "rv", "rv_err"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df[np.isfinite(df["time"]) & np.isfinite(df["rv"]) & np.isfinite(df["rv_err"])].copy()

# Main science set: HR only, split by segment/instrument.
main = df[df["model_inst"].isin(["ESPRESSO_HR_2023_EARLY", "ESPRESSO_HR_2023PLUS"])].copy()

# Center RV per model_inst for periodograms / quick visual QC.
main["rv_centered"] = main["rv"] - main.groupby("model_inst")["rv"].transform("median")

out_full = Path("data_interim/qc/espresso_eso_ccf_curated_hr.csv")
out_pg = OUTDIR / "espresso_eso_ccf_curated_hr_centered.csv"

main.to_csv(out_full, index=False)

main[["time", "rv_centered", "rv_err", "model_inst", "mjd_obs", "berv", "ccf_fwhm", "ccf_bispan", "ccf_contrast"]].rename(
    columns={"rv_centered": "rv"}
).to_csv(out_pg, index=False)

print("Input:", INP)
print("Rows all:", len(df))
print("Rows HR main:", len(main))
print("Saved:", out_full)
print("Saved:", out_pg)

print("\nBy model_inst:")
print(main.groupby("model_inst").agg(
    n=("rv", "size"),
    mjd_min=("mjd_obs", "min"),
    mjd_max=("mjd_obs", "max"),
    rv_median=("rv", "median"),
    rv_std=("rv", "std"),
    rv_centered_std=("rv_centered", "std"),
    err_median=("rv_err", "median"),
).to_string())

print("\nDONE")
