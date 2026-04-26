from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/multinstrument_rvonly_centered.csv")
outpng = Path("data_interim/qc/espresso_only_periodogram_masked.png")

df = pd.read_csv(inp).copy()
df = df[df["instrument"] == "ESPRESSO"].copy()
df = df[df["bjd_tdb"].notna() & df["rv_centered"].notna() & df["rv_err"].notna()].copy()

t = df["bjd_tdb"].values
y = df["rv_centered"].values
e = df["rv_err"].values

min_period = 1.0
max_period = 2000.0
fmin = 1.0 / max_period
fmax = 1.0 / min_period

ls = LombScargle(t, y, e)
freq, power = ls.autopower(
    minimum_frequency=fmin,
    maximum_frequency=fmax,
    samples_per_peak=20
)

period = 1.0 / freq

mask = ~((period > 0.95) & (period < 1.25))
pm = period[mask]
pw = power[mask]

top = np.argsort(pw)[-20:][::-1]

print("\n=== ESPRESSO RV TOP (masked 0.95-1.25 d) ===")
for i in top:
    print(f"P={pm[i]:9.3f} d   power={pw[i]:.6f}")

plt.figure(figsize=(10,6))
plt.plot(period, power, alpha=0.35, label="full")
plt.plot(pm, pw, label="masked 0.95-1.25 d")
for p in [3.24, 11.2, 16.3, 24.99, 36.2, 365.25]:
    plt.axvline(p, linestyle="--", linewidth=1)
plt.xscale("log")
plt.xlabel("Period (days)")
plt.ylabel("Lomb-Scargle power")
plt.title("ESPRESSO-only RV periodogram (daily-alias masked)")
plt.legend()
plt.tight_layout()
plt.savefig(outpng, dpi=160)

print(f"\nSaved {outpng}")
