from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/uves_relative.csv")
out = Path("data_interim/qc/uves_relative_detrended.csv")

df = pd.read_csv(inp).copy()
df = df.sort_values("bjd_tdb").copy()

t = df["bjd_tdb"].values
y = df["y_rel"].values

def local_smooth(t, y, width=15.0):
    ys = np.full_like(y, np.nan, dtype=float)
    for i in range(len(y)):
        dt = np.abs(t - t[i])
        w = np.exp(-0.5 * (dt / width)**2)
        m = np.isfinite(y) & np.isfinite(w)
        if m.sum() < 3:
            continue
        ys[i] = np.sum(w[m] * y[m]) / np.sum(w[m])
    return ys

trend = local_smooth(t, y, width=15.0)
df["y_rel_trend15d"] = trend
df["y_rel_resid15d"] = y - trend

df.to_csv(out, index=False)

print(df[["date_obs","y_rel","y_rel_trend15d","y_rel_resid15d","snr"]].head(25).to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
