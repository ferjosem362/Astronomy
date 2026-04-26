from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/espresso_relative.csv")
out = Path("data_interim/qc/espresso_relative_detrended.csv")

df = pd.read_csv(inp).copy().sort_values("bjd_tdb")

t = df["bjd_tdb"].values
y = df["y_rel"].values

def local_smooth(t, y, width=20.0):
    ys = np.full_like(y, np.nan, dtype=float)
    for i in range(len(y)):
        dt = np.abs(t - t[i])
        w = np.exp(-0.5 * (dt / width)**2)
        m = np.isfinite(y) & np.isfinite(w)
        if m.sum() < 3:
            continue
        ys[i] = np.sum(w[m] * y[m]) / np.sum(w[m])
    return ys

trend = local_smooth(t, y, width=20.0)
df["y_rel_trend20d"] = trend
df["y_rel_resid20d"] = y - trend

df.to_csv(out, index=False)
print(df[["date_obs","y_rel","y_rel_trend20d","y_rel_resid20d"]].head(20).to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
