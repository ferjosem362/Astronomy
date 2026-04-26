from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/multinstrument_rvonly_standardized.csv")
outpng = Path("data_interim/qc/periodogram_compare.png")

df = pd.read_csv(inp).copy()

# tiempos
df = df[df["bjd_tdb"].notna()].copy()

# dataset A: RV física centrada, solo ESPRESSO+HARPS
m_a = df["instrument"].isin(["ESPRESSO","HARPS"])
dfa = df.loc[m_a & df["rv_centered"].notna()].copy()

# dataset B: estandarizado, incluye UVES
m_b = df["instrument"].isin(["ESPRESSO","HARPS","UVES"])
dfb = df.loc[m_b & df["rv_standardized"].notna()].copy()

# arrays
ta = dfa["bjd_tdb"].values
ya = dfa["rv_centered"].values
ea = dfa["rv_err"].fillna(np.nanmedian(dfa["rv_err"])).values

tb = dfb["bjd_tdb"].values
yb = dfb["rv_standardized"].values

# para UVES sus rv_err no son comparables; usamos pesos uniformes en el caso B
# y dejamos que la estandarización mande
min_period = 1.0
max_period = 2000.0

fa_min = 1.0 / max_period
fa_max = 1.0 / min_period

lsa = LombScargle(ta, ya, ea)
freq_a, power_a = lsa.autopower(
    minimum_frequency=fa_min,
    maximum_frequency=fa_max,
    samples_per_peak=15
)

lsb = LombScargle(tb, yb)
freq_b, power_b = lsb.autopower(
    minimum_frequency=fa_min,
    maximum_frequency=fa_max,
    samples_per_peak=15
)

per_a = 1.0 / freq_a
per_b = 1.0 / freq_b

# picos principales
top_a = np.argsort(power_a)[-10:][::-1]
top_b = np.argsort(power_b)[-10:][::-1]

print("\n=== TOP A: ESPRESSO + HARPS ===")
for i in top_a:
    print(f"P={per_a[i]:9.3f} d   power={power_a[i]:.6f}")

print("\n=== TOP B: ESPRESSO + HARPS + UVES(est.) ===")
for i in top_b:
    print(f"P={per_b[i]:9.3f} d   power={power_b[i]:.6f}")

plt.figure(figsize=(10,6))
plt.plot(per_a, power_a, label="ESPRESSO+HARPS")
plt.plot(per_b, power_b, label="ESPRESSO+HARPS+UVES(std)", alpha=0.85)

for p in [3.24, 11.2, 36.2, 365.25]:
    plt.axvline(p, linestyle="--", linewidth=1)

plt.xscale("log")
plt.xlabel("Period (days)")
plt.ylabel("Lomb-Scargle power")
plt.title("Alpha Cen B RV periodogram comparison")
plt.legend()
plt.tight_layout()
plt.savefig(outpng, dpi=160)

print(f"\nSaved {outpng}")
