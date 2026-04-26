from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/espresso_decorrelated.csv")
out = Path("data_interim/qc/espresso_relative.csv")

df = pd.read_csv(inp).copy()
df = df[df["rv_resid_decor"].notna() & df["bjd_tdb"].notna()].copy()

y = df["rv_resid_decor"].astype(float).values
med = np.nanmedian(y)
mad = np.nanmedian(np.abs(y - med))
if not np.isfinite(mad) or mad <= 0:
    mad = np.nanstd(y)
if not np.isfinite(mad) or mad <= 0:
    mad = 1.0

df["y_rel"] = (y - med) / mad
df["y_med"] = med
df["y_scale"] = mad

df.to_csv(out, index=False)
print(df[["date_obs","rv_resid_decor","y_rel","rv_err"]].head(20).to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
