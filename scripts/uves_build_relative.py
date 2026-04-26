from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/uves_pseudorv_sobr_alphacenb_clean_bjd.csv")
out = Path("data_interim/qc/uves_relative.csv")

df = pd.read_csv(inp).copy()
df = df[df["rv_pseudo"].notna() & df["bjd_tdb"].notna()].copy()

y = pd.to_numeric(df["rv_pseudo"], errors="coerce").values
med = np.nanmedian(y)
mad = np.nanmedian(np.abs(y - med))
if not np.isfinite(mad) or mad <= 0:
    mad = np.nanstd(y)
if not np.isfinite(mad) or mad <= 0:
    mad = 1.0

df["instrument"] = "UVES"
df["rv_rel"] = y
df["rv_rel_med"] = med
df["rv_rel_scale"] = mad
df["y_rel"] = (y - med) / mad

df.to_csv(out, index=False)

print(df[["date_obs","rv_pseudo","y_rel","rv_err_pseudo","snr"]].head(20).to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
