from pathlib import Path
import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/espresso_decorrelated.csv")
out = Path("data_interim/qc/espresso_decorrelated_after_2p257d_2p024d.csv")

df = pd.read_csv(inp).copy()
df = df[df["rv_resid_decor"].notna() & df["bjd_tdb"].notna() & df["rv_err"].notna()].copy()

t = df["bjd_tdb"].values
y = df["rv_resid_decor"].values
e = df["rv_err"].values

periods = [2.257, 2.024]
cols = [np.ones_like(t)]
for P0 in periods:
    f0 = 1.0 / P0
    cols.append(np.sin(2*np.pi*f0*t))
    cols.append(np.cos(2*np.pi*f0*t))

X = np.column_stack(cols)
w = 1.0 / (e**2)
sw = np.sqrt(w)

beta, *_ = np.linalg.lstsq(X * sw[:,None], y * sw, rcond=None)

model = X @ beta
resid = y - model

df["model_2p257_2p024"] = model
df["resid_after_2p257_2p024"] = resid
df.to_csv(out, index=False)

ls = LombScargle(t, resid, e)
freq, power = ls.autopower(
    minimum_frequency=1/2000,
    maximum_frequency=1/1.0,
    samples_per_peak=20
)

period = 1.0 / freq
mask = ~((period > 0.95) & (period < 1.25))
pm = period[mask]
pw = power[mask]

top = np.argsort(pw)[-20:][::-1]

print("\n=== TOP after removing 2.257 d and 2.024 d ===")
for i in top:
    print(f"P={pm[i]:9.3f} d   power={pw[i]:.6f}")

print(f"\nSaved {out}")
