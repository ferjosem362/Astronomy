from pathlib import Path
import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/espresso_decorrelated.csv")

df = pd.read_csv(inp).copy()
df = df[df["rv_resid_decor"].notna() & df["bjd_tdb"].notna() & df["rv_err"].notna()].copy()

t = df["bjd_tdb"].values
y = df["rv_resid_decor"].values
e = df["rv_err"].values

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

print("\n=== ESPRESSO DECORRELATED TOP (masked daily) ===")
for i in top:
    print(f"P={pm[i]:9.3f} d   power={pw[i]:.6f}")
