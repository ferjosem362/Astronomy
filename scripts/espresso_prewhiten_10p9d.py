import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

df = pd.read_csv("data_interim/qc/espresso_decorrelated.csv")
df = df[df["rv_resid_decor"].notna() & df["bjd_tdb"].notna() & df["rv_err"].notna()].copy()

t = df["bjd_tdb"].values
y = df["rv_resid_decor"].values
e = df["rv_err"].values

P0 = 10.9
f0 = 1.0 / P0

X = np.column_stack([
    np.ones_like(t),
    np.sin(2*np.pi*f0*t),
    np.cos(2*np.pi*f0*t),
])

w = 1.0 / (e**2)
sw = np.sqrt(w)
beta, *_ = np.linalg.lstsq(X * sw[:,None], y * sw, rcond=None)

resid = y - X @ beta

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

print("\n=== TOP after removing 10.9 d ===")
for i in top:
    print(f"P={pm[i]:9.3f} d   power={pw[i]:.6f}")
