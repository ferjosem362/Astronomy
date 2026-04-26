from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/espresso_nightly_binned.csv")
outpng = Path("data_interim/qc/espresso_nightly_periodogram.png")

df = pd.read_csv(inp).copy()

t = df["bjd_tdb"].values
y = df["rv"].values
e = df["rv_err"].values

fmin = 1.0 / 2000.0
fmax = 1.0 / 1.0

ls = LombScargle(t, y, e)
freq, power = ls.autopower(
    minimum_frequency=fmin,
    maximum_frequency=fmax,
    samples_per_peak=20
)

period = 1.0 / freq
top = np.argsort(power)[-20:][::-1]

print("\n=== ESPRESSO NIGHTLY RV TOP ===")
for i in top:
    print(f"P={period[i]:9.3f} d   power={power[i]:.6f}")

plt.figure(figsize=(10,6))
plt.plot(period, power)
for p in [3.24, 11.2, 12.48, 16.3, 24.99, 36.2, 53.3, 217.0, 365.25]:
    plt.axvline(p, linestyle="--", linewidth=1)
plt.xscale("log")
plt.xlabel("Period (days)")
plt.ylabel("Lomb-Scargle power")
plt.title("ESPRESSO nightly-binned RV periodogram")
plt.tight_layout()
plt.savefig(outpng, dpi=160)

print(f"\nSaved {outpng}")
