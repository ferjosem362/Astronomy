from pathlib import Path
from astropy.io import fits
import pandas as pd

root = Path("data_raw/harps")
rows = []

for fp in sorted(root.glob("*.fits")):
    try:
        hdul = fits.open(fp)
        hdr0 = hdul[0].header
        hdr1 = hdul[1].header if len(hdul) > 1 else {}

        rows.append({
            "file": str(fp),
            "filename": fp.name,
            "object": hdr0.get("OBJECT"),
            "instrume": hdr0.get("INSTRUME"),
            "date_obs": hdr0.get("DATE-OBS"),
            "mjd_obs": hdr0.get("MJD-OBS"),
            "exptime": hdr0.get("EXPTIME"),
            "snr": hdr0.get("SNR") or hdr0.get("HIERARCH ESO DRS SPE EXT SN50") or hdr1.get("SNR"),
            "berv": hdr0.get("HIERARCH ESO DRS BERV") or hdr1.get("HIERARCH ESO DRS BERV"),
            "rv": hdr0.get("HIERARCH ESO DRS CCF RV") or hdr1.get("HIERARCH ESO DRS CCF RV"),
            "rv_err": hdr0.get("HIERARCH ESO DRS CCF NOISE") or hdr1.get("HIERARCH ESO DRS CCF NOISE"),
            "bis": hdr0.get("HIERARCH ESO DRS BIS SPAN") or hdr1.get("HIERARCH ESO DRS BIS SPAN"),
            "fwhm": hdr0.get("HIERARCH ESO DRS CCF FWHM") or hdr1.get("HIERARCH ESO DRS CCF FWHM"),
            "airmass": hdr0.get("HIERARCH ESO TEL AIRM START"),
            "prog_id": hdr0.get("HIERARCH ESO OBS PROG ID"),
            "arcfile": hdr0.get("ARCFILE"),
            "origfile": hdr0.get("ORIGFILE"),
            "n_hdu": len(hdul),
        })
        hdul.close()
    except Exception as e:
        rows.append({"file": str(fp), "error": str(e)})

df = pd.DataFrame(rows)
out = Path("data_interim/qc/harps_master_minimal.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
print(df.to_string(index=False))
print(f"\nSaved: {out}")
