from pathlib import Path
from astropy.io import fits
import pandas as pd

rows = []

for ccf_fp in sorted(Path("data_interim/harps_aux_unpacked").rglob("*_ccf_G2_A.fits")):
    stem = ccf_fp.name.replace("_ccf_G2_A.fits", "")
    folder = ccf_fp.parent

    bis_fp = folder / f"{stem}_bis_G2_A.fits"
    s1d_fp = folder / f"{stem}_s1d_A.fits"
    e2ds_fp = folder / f"{stem}_e2ds_A.fits"

    try:
        with fits.open(ccf_fp) as hdul:
            hdr = hdul[0].header

            rows.append({
                "stem": stem,
                "ccf_file": str(ccf_fp),
                "bis_file": str(bis_fp) if bis_fp.exists() else None,
                "s1d_file": str(s1d_fp) if s1d_fp.exists() else None,
                "e2ds_file": str(e2ds_fp) if e2ds_fp.exists() else None,
                "object": hdr.get("OBJECT"),
                "date_obs": hdr.get("DATE-OBS"),
                "mjd_obs": hdr.get("MJD-OBS"),
                "bjd": hdr.get("ESO DRS BJD"),
                "berv": hdr.get("ESO DRS BERV"),
                "rv": hdr.get("ESO DRS CCF RV"),
                "rv_err": hdr.get("ESO DRS CCF NOISE"),
                "fwhm": hdr.get("ESO DRS CCF FWHM"),
                "contrast": hdr.get("ESO DRS CCF CONTRAST"),
                "mask": hdr.get("ESO DRS CCF MASK"),
                "airmass": hdr.get("HIERARCH ESO TEL AIRM START"),
                "snr": hdr.get("HIERARCH ESO DRS SPE EXT SN50") or hdr.get("SNR"),
                "prog_id": hdr.get("HIERARCH ESO OBS PROG ID"),
            })
    except Exception as e:
        rows.append({"ccf_file": str(ccf_fp), "error": str(e)})

df = pd.DataFrame(rows).sort_values("date_obs")
out = Path("data_interim/qc/harps_aux_master.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print(df.to_string(index=False))
print(f"\nSaved: {out}")
