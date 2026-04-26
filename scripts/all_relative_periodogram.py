import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

df = pd.read_csv("data_interim/qc/all_relative_clean.csv")
df = df[df["y_clean"].notna() & df["bjd_tdb"].notna()].copy()

t = df["bjd_tdb"].values
y = df["y_clean"].values
e = np.ones_like(y)

ls = LombScargle(t, y, e)
freq, power = ls.autopower(
    minimum_frequency=1/3000,
    maximum_frequency=1/1.0,
    samples_per_peak=20
)

period = 1.0 / freq
mask = ~((period > 0.95) & (period < 1.25))
pm = period[mask]
pw = power[mask]

top = np.argsort(pw)[-30:][::-1]

print("\n=== ALL RELATIVE CLEAN TOP ===")
for i in top:
    print(f"P={pm[i]:10.3f} d   power={pw[i]:.6f}")
