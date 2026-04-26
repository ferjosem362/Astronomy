from pathlib import Path
import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/espresso_relative_detrended.csv")
df = pd.read_csv(inp).copy()
df = df[df["y_rel_resid20d"].notna() & df["bjd_tdb"].notna()].copy()

t = df["bjd_tdb"].values
y = df["y_rel_resid20d"].values
e = pd.to_numeric(df["rv_err"], errors="coerce").values
bad = ~np.isfinite(e) | (e <= 0)
if bad.any():
    e[bad] = np.nanmedian(e[~bad]) if (~bad).any() else 1.0

ls = LombScargle(t, y, e)
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

print("\n=== ESPRESSO RELATIVE DETRENDED TOP ===")
for i in top:
    print(f"P={pm[i]:9.3f} d   power={pw[i]:.6f}")
