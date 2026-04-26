from pathlib import Path
from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

files = sorted(Path("data_raw/harps").glob("*.fits"))

specs = []
for fp in files:
    with fits.open(fp) as hdul:
        hdr = hdul[0].header
        tab = hdul[1].data
        wave = np.array(tab["WAVE"][0], dtype=float)
        flux = np.array(tab["FLUX"][0], dtype=float)

        # ventana pequeña y limpia para prueba
        m = (wave > 5000) & (wave < 5100) & np.isfinite(flux)
        specs.append({
            "file": fp.name,
            "bjd": hdr.get("ESO DRS BJD"),
            "berv": hdr.get("ESO DRS BERV"),
            "wave": wave[m],
            "flux": flux[m],
        })

# malla común
wmin = max(s["wave"].min() for s in specs)
wmax = min(s["wave"].max() for s in specs)
grid = np.linspace(wmin, wmax, 5000)

F = []
meta = []
for s in specs:
    y = np.interp(grid, s["wave"], s["flux"])
    # normalización tosca
    y = y / np.nanmedian(y)
    F.append(y)
    meta.append((s["file"], s["bjd"], s["berv"]))

F = np.array(F)
template = np.nanmedian(F, axis=0)
resid = F - template

# RMS por exposición
rms = np.sqrt(np.nanmean(resid**2, axis=1))

df = pd.DataFrame(meta, columns=["file","bjd","berv"])
df["resid_rms"] = rms
print(df.to_string(index=False))

plt.figure(figsize=(8,4))
for i in range(min(5, len(F))):
    plt.plot(grid, F[i], alpha=0.7, label=meta[i][0][:22])
plt.xlabel("Wavelength [A]")
plt.ylabel("Normalized flux")
plt.title("Mini-LBL window 5000-5100 A")
plt.tight_layout()
plt.show()

plt.figure(figsize=(6,4))
plt.scatter(df["berv"], df["resid_rms"])
plt.xlabel("BERV")
plt.ylabel("Residual RMS")
plt.title("Residual RMS vs BERV")
plt.tight_layout()
plt.show()
