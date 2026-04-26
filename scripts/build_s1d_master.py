from pathlib import Path
from astropy.io import fits
import pandas as pd
import numpy as np

rows = []

for fp in sorted(Path("data_raw/harps").glob("*.fits")):
    try:
        with fits.open(fp) as hdul:
            hdr = hdul[0].header
            tab = hdul[1].data

            wave = np.array(tab["WAVE"][0], dtype=float)
            flux = np.array(tab["FLUX"][0], dtype=float)
            err  = np.array(tab["ERR"][0], dtype=float)

            finite_flux = np.isfinite(flux)
            finite_err  = np.isfinite(err)

            rows.append({
                "file": str(fp),
                "stem": fp.name.replace("_s1d_A_DRS_HARPS_3.5_ESOSDP.fits",""),
                "object": hdr.get("OBJECT"),
                "date_obs": hdr.get("DATE-OBS"),
                "mjd_obs": hdr.get("MJD-OBS"),
                "bjd": hdr.get("ESO DRS BJD"),
                "berv": hdr.get("ESO DRS BERV"),
                "snr": hdr.get("HIERARCH ESO DRS SPE EXT SN50") or hdr.get("SNR"),
                "airmass": hdr.get("HIERARCH ESO TEL AIRM START"),
                "wave_min": float(np.nanmin(wave)),
                "wave_max": float(np.nanmax(wave)),
                "n_pix": int(len(wave)),
                "finite_flux_frac": float(finite_flux.mean()),
                "finite_err_frac": float(finite_err.mean()),
                "median_flux": float(np.nanmedian(flux)),
                "median_err": float(np.nanmedian(err)) if finite_err.any() else None,
            })
    except Exception as e:
        rows.append({"file": str(fp), "error": str(e)})

df = pd.DataFrame(rows).sort_values("bjd")
out = Path("data_interim/qc/s1d_master.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
print(df.to_string(index=False))
print(f"\nSaved: {out}")
