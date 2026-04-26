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
            hccf = hdul[0].header

        bis_rv = None
        bis_span = None
        if bis_fp.exists():
            with fits.open(bis_fp) as hdul:
                hbis = hdul[0].header
                bis_rv = hbis.get("ESO DRS BIS RV")
                bis_span = hbis.get("ESO DRS BIS SPAN")

        rows.append({
            "stem": stem,
            "object": hccf.get("OBJECT"),
            "date_obs": hccf.get("DATE-OBS"),
            "mjd_obs": hccf.get("MJD-OBS"),
            "bjd": hccf.get("ESO DRS BJD"),
            "berv": hccf.get("ESO DRS BERV"),
            "rv": hccf.get("ESO DRS CCF RV"),
            "rv_err": hccf.get("ESO DRS CCF NOISE"),
            "fwhm": hccf.get("ESO DRS CCF FWHM"),
            "contrast": hccf.get("ESO DRS CCF CONTRAST"),
            "mask": hccf.get("ESO DRS CCF MASK"),
            "bis_rv": bis_rv,
            "bis_span": bis_span,
            "snr": hccf.get("HIERARCH ESO DRS SPE EXT SN50") or hccf.get("SNR"),
            "airmass": hccf.get("HIERARCH ESO TEL AIRM START"),
            "prog_id": hccf.get("HIERARCH ESO OBS PROG ID"),
            "ccf_file": str(ccf_fp),
            "bis_file": str(bis_fp) if bis_fp.exists() else None,
            "s1d_file": str(s1d_fp) if s1d_fp.exists() else None,
            "e2ds_file": str(e2ds_fp) if e2ds_fp.exists() else None,
        })
    except Exception as e:
        rows.append({"ccf_file": str(ccf_fp), "error": str(e)})

df = pd.DataFrame(rows).sort_values("bjd")
out = Path("data_interim/qc/harps_aux_master_with_bis.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
print(df.to_string(index=False))
print(f"\nSaved: {out}")
