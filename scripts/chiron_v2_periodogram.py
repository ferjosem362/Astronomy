from pathlib import Path
import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/chiron_v2_decorrelated.csv")
df = pd.read_csv(inp).copy()
df = df[df["rv_resid_decor"].notna() & df["bjd_tdb"].notna()].copy()

t = df["bjd_tdb"].values
y = df["rv_resid_decor"].values
e = pd.to_numeric(df["rv_rel_err"], errors="coerce").values
bad = ~np.isfinite(e) | (e <= 0)
if bad.any():
    e[bad] = np.nanmedian(e[~bad]) if (~bad).any() else 1.0

ls = LombScargle(t, y, e)
freq, power = ls.autopower(
    minimum_frequency=1/1000,
    maximum_frequency=1/1.0,
    samples_per_peak=20
)

period = 1.0 / freq
mask = ~((period > 0.95) & (period < 1.25))
pm = period[mask]
pw = power[mask]

top = np.argsort(pw)[-20:][::-1]

print("\n=== CHIRON v2 DECORRELATED TOP ===")
for i in top:
    print(f"P={pm[i]:9.3f} d   power={pw[i]:.6f}")
