from pathlib import Path
from astropy.io import fits
import pandas as pd
import numpy as np
import re

base = Path("data_raw/espresso_sample10")
files = sorted(base.glob("*.fits"))

def classify(name: str) -> str:
    if "SCCA" in name:
        return "ccf"
    if "SDMB" in name:
        return "drift"
    if "SFLX" in name or "SOBF" in name:
        return "spectrum"
    return "other"

def obskey(name: str) -> str:
    # quita prefijo de tipo de producto y extensión
    x = name.replace(".fits", "")
    x = re.sub(r"^ES_[A-Z]+_", "", x)
    return x

groups = {}
for fp in files:
    key = obskey(fp.name)
    groups.setdefault(key, []).append(fp)

rows = []

for key, fps in groups.items():
    row = {
        "instrument": "ESPRESSO",
        "stem": key,
        "target_name": None,
        "date_obs": None,
        "bjd_tdb": np.nan,
        "mjd_obs": np.nan,
        "rv": np.nan,
        "rv_err": np.nan,
        "berv": np.nan,
        "fwhm": np.nan,
        "contrast": np.nan,
        "bis_span": np.nan,
        "snr": np.nan,
        "airmass": np.nan,
        "drift": np.nan,
        "aca_sep_arcsec": np.nan,
        "aca_contam_proxy": np.nan,
        "file_spectrum": None,
        "file_ccf": None,
        "file_bis": None,
        "file_e2ds": None,
        "file_drift": None,
        "quality_flag": 0,
        "notes": "",
    }

    for fp in fps:
        kind = classify(fp.name)
        try:
            with fits.open(fp) as hdul:
                hdr = hdul[0].header

                if row["target_name"] is None:
                    row["target_name"] = hdr.get("OBJECT")
                if row["date_obs"] is None:
                    row["date_obs"] = hdr.get("DATE-OBS")
                if pd.isna(row["mjd_obs"]):
                    row["mjd_obs"] = hdr.get("MJD-OBS", np.nan)

                for k in ["ESO QC BJD", "ESO DRS BJD", "BJD"]:
                    if k in hdr and pd.isna(row["bjd_tdb"]):
                        row["bjd_tdb"] = hdr[k]

                for k in ["ESO QC BERV", "ESO DRS BERV", "BERV"]:
                    if k in hdr and pd.isna(row["berv"]):
                        row["berv"] = hdr[k]

                if pd.isna(row["airmass"]):
                    row["airmass"] = hdr.get("HIERARCH ESO TEL AIRM START", np.nan)

                # SNR útil
                for k in ["ESO QC ORDER110 SNR", "ESO QC ORDER55 SNR", "SNR"]:
                    if k in hdr and pd.isna(row["snr"]):
                        row["snr"] = hdr[k]

                if kind == "ccf":
                    row["file_ccf"] = str(fp)

                    for k in [
                        "ESO QC CCF RV",
                        "HIERARCH ESO QC CCF RV",
                        "ESO DRS CCF RV",
                        "RV",
                    ]:
                        if k in hdr and pd.isna(row["rv"]):
                            row["rv"] = hdr[k]

                    for k in [
                        "ESO QC CCF RV ERROR",
                        "HIERARCH ESO QC CCF RV ERROR",
                        "ESO QC CCF RV ERR",
                        "ESO DRS CCF NOISE",
                        "RVERR",
                    ]:
                        if k in hdr and pd.isna(row["rv_err"]):
                            row["rv_err"] = hdr[k]

                    for k in ["ESO QC CCF FWHM", "HIERARCH ESO QC CCF FWHM", "ESO DRS CCF FWHM", "FWHM"]:
                        if k in hdr and pd.isna(row["fwhm"]):
                            row["fwhm"] = hdr[k]

                    for k in ["ESO QC CCF CONTRAST", "HIERARCH ESO QC CCF CONTRAST", "ESO DRS CCF CONTRAST", "CONTRAST"]:
                        if k in hdr and pd.isna(row["contrast"]):
                            row["contrast"] = hdr[k]

                    for k in ["ESO QC BIS SPAN", "HIERARCH ESO QC BIS SPAN", "ESO DRS BIS SPAN", "BISSPAN"]:
                        if k in hdr and pd.isna(row["bis_span"]):
                            row["bis_span"] = hdr[k]

                elif kind == "drift":
                    row["file_drift"] = str(fp)

                    # dos detectores: media simple
                    d0 = hdr.get("ESO QC DRIFT DET0 MEAN", np.nan)
                    d1 = hdr.get("ESO QC DRIFT DET1 MEAN", np.nan)
                    vals = [x for x in [d0, d1] if pd.notna(x)]
                    if vals and pd.isna(row["drift"]):
                        row["drift"] = float(np.mean(vals))

                elif kind == "spectrum":
                    if row["file_spectrum"] is None:
                        row["file_spectrum"] = str(fp)

        except Exception as e:
            row["notes"] += f"{fp.name}:{e}; "

    rows.append(row)

df = pd.DataFrame(rows).sort_values("bjd_tdb")
out = Path("data_interim/qc/espresso_master.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print(f"saved {out} rows={len(df)}")
print(df[[
    "instrument","stem","target_name","date_obs","bjd_tdb",
    "rv","rv_err","berv","fwhm","contrast","bis_span","snr","drift"
]].to_string(index=False))
