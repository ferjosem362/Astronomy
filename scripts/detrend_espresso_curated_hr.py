#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd

INP = Path("data_interim/qc/espresso_eso_ccf_curated_hr.csv")
OUT1 = Path("data_interim/qc/espresso_eso_ccf_curated_hr_detrended.csv")
OUT2 = Path("data_interim/qc/periodogram_inputs/espresso_eso_ccf_curated_hr_detrended.csv")
LOG = Path("logs/detrend_espresso_curated_hr_summary.txt")

df = pd.read_csv(INP)

for c in ["time", "mjd_obs", "rv", "rv_err", "berv", "ccf_fwhm", "ccf_bispan", "ccf_contrast"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

df = df[np.isfinite(df["time"]) & np.isfinite(df["rv"]) & np.isfinite(df["rv_err"])].copy()
df = df.sort_values(["model_inst", "time"]).reset_index(drop=True)

rows = []
parts = []

for inst, g in df.groupby("model_inst", sort=False):
    g = g.copy()
    t = g["time"].to_numpy()
    y = g["rv"].to_numpy()
    e = g["rv_err"].to_numpy()

    # Para no comernos señales de días, quitamos solo tendencia larga.
    # Grado 2 si hay suficiente baseline; grado 1 si no.
    baseline = np.nanmax(t) - np.nanmin(t)
    deg = 2 if (len(g) >= 50 and baseline > 100) else 1

    x = t - np.nanmedian(t)
    w = 1.0 / np.maximum(e, np.nanmedian(e))**2

    coef = np.polyfit(x, y, deg=deg, w=np.sqrt(w))
    trend = np.polyval(coef, x)

    g["rv_trend"] = trend
    g["rv_detrended"] = y - trend
    g["rv_detrended_centered"] = g["rv_detrended"] - np.nanmedian(g["rv_detrended"])

    rows.append({
        "model_inst": inst,
        "n": len(g),
        "bjd_min": float(np.nanmin(t)),
        "bjd_max": float(np.nanmax(t)),
        "baseline_d": float(baseline),
        "poly_degree": deg,
        "rv_raw_std": float(np.nanstd(y, ddof=1)),
        "rv_detrended_std": float(np.nanstd(g["rv_detrended_centered"], ddof=1)),
        "rv_err_median": float(np.nanmedian(e)),
        "coef": repr(coef.tolist()),
    })
    parts.append(g)

out = pd.concat(parts, ignore_index=True)
summary = pd.DataFrame(rows)

OUT1.parent.mkdir(parents=True, exist_ok=True)
OUT2.parent.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)

out.to_csv(OUT1, index=False)

out[[
    "time",
    "rv_detrended_centered",
    "rv_err",
    "model_inst",
    "mjd_obs",
    "berv",
    "ccf_fwhm",
    "ccf_bispan",
    "ccf_contrast",
]].rename(columns={"rv_detrended_centered": "rv"}).to_csv(OUT2, index=False)

with open(LOG, "w") as f:
    f.write(summary.to_string(index=False))
    f.write("\n")

print("Input:", INP)
print("Rows:", len(out))
print("Saved:", OUT1)
print("Saved:", OUT2)
print("Saved:", LOG)
print()
print(summary.to_string(index=False))
print("\nDONE")
