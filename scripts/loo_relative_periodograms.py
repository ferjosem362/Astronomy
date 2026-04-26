import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

df = pd.read_csv("data_interim/qc/all_relative_clean.csv")
df = df[df["y_clean"].notna() & df["bjd_tdb"].notna()].copy()

instruments = sorted(df["instrument"].dropna().unique())

for leave_out in instruments:
    d = df[df["instrument"] != leave_out].copy()
    t = d["bjd_tdb"].values
    y = d["y_clean"].values
    e = np.ones_like(y)

    ls = LombScargle(t, y, e)
    freq, power = ls.autopower(
        minimum_frequency=1/5000,
        maximum_frequency=1/1.0,
        samples_per_peak=20
    )

    period = 1.0 / freq
    mask = ~((period > 0.95) & (period < 1.25))
    pm = period[mask]
    pw = power[mask]
    top = np.argsort(pw)[-15:][::-1]

    print(f"\n=== LEAVE OUT {leave_out} ===")
    for i in top:
        print(f"P={pm[i]:10.3f} d   power={pw[i]:.6f}")
