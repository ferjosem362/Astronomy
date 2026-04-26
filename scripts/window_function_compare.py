from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/multinstrument_rvonly_standardized.csv")
outpng = Path("data_interim/qc/window_function_compare.png")

df = pd.read_csv(inp).copy()
df = df[df["bjd_tdb"].notna()].copy()

# A: ESPRESSO + HARPS
dfa = df[df["instrument"].isin(["ESPRESSO","HARPS"])].copy()
ta = dfa["bjd_tdb"].values
ya = np.ones_like(ta)

# B: ESPRESSO + HARPS + UVES
dfb = df[df["instrument"].isin(["ESPRESSO","HARPS","UVES"])].copy()
tb = dfb["bjd_tdb"].values
yb = np.ones_like(tb)

min_period = 1.0
max_period = 2000.0
fmin = 1.0 / max_period
fmax = 1.0 / min_period

lsa = LombScargle(ta, ya, center_data=True)
freq_a, power_a = lsa.autopower(
    minimum_frequency=fmin,
    maximum_frequency=fmax,
    samples_per_peak=15
)

lsb = LombScargle(tb, yb, center_data=True)
freq_b, power_b = lsb.autopower(
    minimum_frequency=fmin,
    maximum_frequency=fmax,
    samples_per_peak=15
)

per_a = 1.0 / freq_a
per_b = 1.0 / freq_b

top_a = np.argsort(power_a)[-10:][::-1]
top_b = np.argsort(power_b)[-10:][::-1]

print("\n=== WINDOW TOP A: ESPRESSO + HARPS ===")
for i in top_a:
    print(f"P={per_a[i]:9.3f} d   power={power_a[i]:.6f}")

print("\n=== WINDOW TOP B: ESPRESSO + HARPS + UVES ===")
for i in top_b:
    print(f"P={per_b[i]:9.3f} d   power={power_b[i]:.6f}")

plt.figure(figsize=(10,6))
plt.plot(per_a, power_a, label="Window: ESPRESSO+HARPS")
plt.plot(per_b, power_b, label="Window: ESPRESSO+HARPS+UVES", alpha=0.85)

for p in [1.0, 1.084, 3.24, 12.45, 15.2, 24.99, 36.2, 365.25]:
    plt.axvline(p, linestyle="--", linewidth=1)

plt.xscale("log")
plt.xlabel("Period (days)")
plt.ylabel("Window power")
plt.title("Sampling window function")
plt.legend()
plt.tight_layout()
plt.savefig(outpng, dpi=160)

print(f"\nSaved {outpng}")
