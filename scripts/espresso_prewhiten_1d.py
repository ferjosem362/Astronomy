from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

inp = Path("data_interim/qc/espresso_nightly_binned.csv")
outcsv = Path("data_interim/qc/espresso_nightly_residuals_after_1p106d.csv")
outpng = Path("data_interim/qc/espresso_nightly_periodogram_after_1p106d.png")

df = pd.read_csv(inp).copy()

t = df["bjd_tdb"].values
y = df["rv"].values
e = df["rv_err"].values

# periodo a sustraer
P0 = 1.106
f0 = 1.0 / P0

# diseño lineal: offset + sin + cos
w = 1.0 / (e**2)
X = np.column_stack([
    np.ones_like(t),
    np.sin(2*np.pi*f0*t),
    np.cos(2*np.pi*f0*t),
])

W = np.diag(w)
beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y)
model = X @ beta
resid = y - model

df["rv_model_1p106d"] = model
df["rv_resid"] = resid
df.to_csv(outcsv, index=False)

# periodograma de residuos
ls = LombScargle(t, resid, e)
freq, power = ls.autopower(
    minimum_frequency=1.0/2000.0,
    maximum_frequency=1.0/1.0,
    samples_per_peak=20
)

period = 1.0 / freq
top = np.argsort(power)[-20:][::-1]

print("\n=== RESIDUAL TOP after removing 1.106 d ===")
for i in top:
    print(f"P={period[i]:9.3f} d   power={power[i]:.6f}")

plt.figure(figsize=(10,6))
plt.plot(period, power)
for p in [3.24, 11.2, 12.48, 16.3, 24.99, 36.2, 53.3, 217.0, 365.25]:
    plt.axvline(p, linestyle="--", linewidth=1)
plt.xscale("log")
plt.xlabel("Period (days)")
plt.ylabel("Lomb-Scargle power")
plt.title("ESPRESSO nightly periodogram after removing 1.106 d")
plt.tight_layout()
plt.savefig(outpng, dpi=160)

print(f"\nSaved {outcsv}")
print(f"Saved {outpng}")
