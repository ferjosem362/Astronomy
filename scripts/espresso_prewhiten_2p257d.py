from pathlib import Path
import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/espresso_decorrelated.csv")
out = Path("data_interim/qc/espresso_decorrelated_after_2p257d.csv")

df = pd.read_csv(inp).copy()
df = df[df["rv_resid_decor"].notna() & df["bjd_tdb"].notna() & df["rv_err"].notna()].copy()

t = df["bjd_tdb"].values
y = df["rv_resid_decor"].values
e = df["rv_err"].values

P0 = 2.257
f0 = 1.0 / P0

w = 1.0 / (e**2)
X = np.column_stack([
    np.ones_like(t),
    np.sin(2*np.pi*f0*t),
    np.cos(2*np.pi*f0*t),
])
sw = np.sqrt(w)
beta, *_ = np.linalg.lstsq(X * sw[:,None], y * sw, rcond=None)

model = X @ beta
resid = y - model

df["model_2p257d"] = model
df["resid_after_2p257d"] = resid
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

print("\n=== TOP after removing 2.257 d ===")
for i in top:
    print(f"P={pm[i]:9.3f} d   power={pw[i]:.6f}")

print(f"\nSaved {out}")
