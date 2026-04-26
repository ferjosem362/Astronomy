from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/multinstrument_rvonly_centered.csv")
outpng = Path("data_interim/qc/espresso_only_window.png")

df = pd.read_csv(inp).copy()
df = df[df["instrument"] == "ESPRESSO"].copy()
df = df[df["bjd_tdb"].notna()].copy()

t = df["bjd_tdb"].values
y = np.ones(len(t), dtype=float)

min_period = 1.0
max_period = 2000.0
fmin = 1.0 / max_period
fmax = 1.0 / min_period

ls = LombScargle(t, y, fit_mean=False, center_data=False, normalization="psd")
freq, power = ls.autopower(
    minimum_frequency=fmin,
    maximum_frequency=fmax,
    samples_per_peak=20
)

period = 1.0 / freq

good = np.isfinite(power)
period = period[good]
power = power[good]

top = np.argsort(power)[-15:][::-1]

print("\n=== ESPRESSO WINDOW TOP ===")
for i in top:
    print(f"P={period[i]:9.3f} d   power={power[i]:.6f}")

plt.figure(figsize=(10,6))
plt.plot(period, power)
for p in [1.0, 3.24, 11.2, 16.3, 24.99, 36.2, 365.25]:
    plt.axvline(p, linestyle="--", linewidth=1)
plt.xscale("log")
plt.xlabel("Period (days)")
plt.ylabel("Window power")
plt.title("ESPRESSO-only sampling window")
plt.tight_layout()
plt.savefig(outpng, dpi=160)

print(f"\nSaved {outpng}")
