from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/uves_relative_detrended.csv")
out = Path("data_interim/qc/uves_relative_decorrelated.csv")

df = pd.read_csv(inp).copy()
df = df[df["y_rel_resid15d"].notna()].copy()

Xcols = []
if "snr" in df.columns:
    snr = pd.to_numeric(df["snr"], errors="coerce")
    if np.isfinite(np.nanstd(snr)) and np.nanstd(snr) > 0:
        df["snr_z"] = (snr - np.nanmean(snr)) / np.nanstd(snr)
        Xcols.append("snr_z")

y = df["y_rel_resid15d"].astype(float).values

if Xcols:
    A = np.column_stack([np.ones(len(df))] + [df[c].values for c in Xcols])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    model = A @ beta
else:
    beta = np.array([np.nan])
    model = np.nanmean(y) * np.ones_like(y)

df["y_rel_model"] = model
df["y_rel_decor"] = y - model

df.to_csv(out, index=False)

print("Used features:", Xcols)
print("Beta:", beta)
print(df[["date_obs","y_rel_resid15d","y_rel_model","y_rel_decor","snr"]].head(25).to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
