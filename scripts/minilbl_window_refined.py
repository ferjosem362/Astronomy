from pathlib import Path
from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

W1, W2 = 5035.0, 5055.0

files = sorted(Path("data_raw/harps").glob("*.fits"))
specs = []

for fp in files:
    with fits.open(fp) as hdul:
        hdr = hdul[0].header
        tab = hdul[1].data
        wave = np.array(tab["WAVE"][0], dtype=float)
        flux = np.array(tab["FLUX"][0], dtype=float)

    m = (wave > W1) & (wave < W2) & np.isfinite(flux)
    w = wave[m]
    f = flux[m]

    # pseudo-continuo lineal con percentiles altos
    p = np.isfinite(f)
    x = w[p]
    y = f[p]
    thresh = np.nanpercentile(y, 80)
    q = y >= thresh

    if q.sum() >= 2:
        coeff = np.polyfit(x[q], y[q], 1)
        cont = np.polyval(coeff, w)
        fn = f / cont
    else:
        fn = f / np.nanmedian(f)

    # recorte suave de extremos
    lo, hi = np.nanpercentile(fn, [1, 99])
    fn = np.clip(fn, lo, hi)

    specs.append({
        "file": fp.name,
        "bjd": hdr.get("ESO DRS BJD"),
        "berv": hdr.get("ESO DRS BERV"),
        "wave": w,
        "flux": fn,
    })

wmin = max(s["wave"].min() for s in specs)
wmax = min(s["wave"].max() for s in specs)
grid = np.linspace(wmin, wmax, 2500)

F = []
meta = []
for s in specs:
    y = np.interp(grid, s["wave"], s["flux"])
    F.append(y)
    meta.append((s["file"], s["bjd"], s["berv"]))

F = np.array(F)
template = np.nanmedian(F, axis=0)
resid = F - template
rms = np.sqrt(np.nanmean(resid**2, axis=1))

df = pd.DataFrame(meta, columns=["file","bjd","berv"])
df["resid_rms"] = rms
print(df.to_string(index=False))

plt.figure(figsize=(8,4))
for i in range(min(5, len(F))):
    plt.plot(grid, F[i], alpha=0.7)
plt.xlabel("Wavelength [A]")
plt.ylabel("Normalized flux")
plt.title(f"Refined mini-LBL {W1:.1f}-{W2:.1f} A")
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,4))
plt.scatter(df["berv"], df["resid_rms"])
plt.xlabel("BERV")
plt.ylabel("Residual RMS")
plt.title("Residual RMS vs BERV")
plt.tight_layout()
plt.show()
