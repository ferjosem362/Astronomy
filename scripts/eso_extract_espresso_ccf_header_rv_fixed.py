#!/usr/bin/env python3
from pathlib import Path
import re
import numpy as np
import pandas as pd
from astropy.io import fits

INDIR = Path("data_raw/espresso/ccf")
OUTDIR = Path("data_processed/espresso")
METADIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)
METADIR.mkdir(parents=True, exist_ok=True)

files = sorted(INDIR.glob("*.fits"))
if not files:
    files = sorted(Path("data_raw/espresso").glob("**/*SCCA*.fits"))

def get_header_value(header, names, default=np.nan):
    for name in names:
        if name in header:
            return header.get(name)
    # astropy sometimes normalizes HIERARCH cards to bare ESO...
    upper_map = {str(k).upper(): k for k in header.keys()}
    for name in names:
        k = upper_map.get(str(name).upper())
        if k is not None:
            return header.get(k)
    return default

def to_float(x):
    try:
        if x is None:
            return np.nan
        return float(x)
    except Exception:
        return np.nan

def parse_filename(fp):
    s = fp.name
    # ES_SCCA_3548074_2023-04-08T05_28_09.598_HR_1x1_U1__ADP....
    m = re.search(r"ES_SCCA_(\d+)_(\d{4}-\d{2}-\d{2}T\d{2}_\d{2}_\d{2}\.\d+)_([^_]+)_", s)
    if not m:
        return "", "", ""
    obs_id = m.group(1)
    date_from_name = m.group(2).replace("_", ":", 2)
    # fix only time underscores
    date_from_name = m.group(2)
    date_from_name = date_from_name[:13] + ":" + date_from_name[14:16] + ":" + date_from_name[17:]
    mode_from_name = m.group(3)
    return obs_id, date_from_name, mode_from_name

rows = []

for fp in files:
    row = {
        "filename": fp.name,
        "file": str(fp),
        "ok": False,
        "error": "",
    }
    try:
        with fits.open(fp, memmap=True) as hdul:
            h = hdul[0].header

            obs_id, date_from_name, mode_from_name = parse_filename(fp)

            date_obs = get_header_value(h, ["DATE-OBS"], "")
            mjd_obs = to_float(get_header_value(h, ["MJD-OBS"]))
            exptime = to_float(get_header_value(h, ["EXPTIME"]))
            obj = get_header_value(h, ["OBJECT", "ESO OBS TARG NAME"], "")
            target = get_header_value(h, ["ESO OBS TARG NAME", "HIERARCH ESO OBS TARG NAME"], "")
            ins_mode = get_header_value(h, ["ESO INS MODE", "HIERARCH ESO INS MODE"], "")
            pro_catg = get_header_value(h, ["ESO PRO CATG", "HIERARCH ESO PRO CATG"], "")
            prodcatg = get_header_value(h, ["PRODCATG"], "")

            rv = to_float(get_header_value(h, ["ESO QC CCF RV", "HIERARCH ESO QC CCF RV"]))
            rv_err = to_float(get_header_value(h, ["ESO QC CCF RV ERROR", "HIERARCH ESO QC CCF RV ERROR"]))
            berv = to_float(get_header_value(h, ["ESO QC BERV", "HIERARCH ESO QC BERV"]))
            fwhm = to_float(get_header_value(h, ["ESO QC CCF FWHM", "HIERARCH ESO QC CCF FWHM"]))
            fwhm_err = to_float(get_header_value(h, ["ESO QC CCF FWHM ERROR", "HIERARCH ESO QC CCF FWHM ERROR"]))
            contrast = to_float(get_header_value(h, ["ESO QC CCF CONTRAST", "HIERARCH ESO QC CCF CONTRAST"]))
            contrast_err = to_float(get_header_value(h, ["ESO QC CCF CONTRAST ERROR", "HIERARCH ESO QC CCF CONTRAST ERROR"]))
            bis = to_float(get_header_value(h, ["ESO QC CCF BIS SPAN", "HIERARCH ESO QC CCF BIS SPAN"]))
            bis_err = to_float(get_header_value(h, ["ESO QC CCF BIS SPAN ERROR", "HIERARCH ESO QC CCF BIS SPAN ERROR"]))
            mask = get_header_value(h, ["ESO QC CCF MASK", "HIERARCH ESO QC CCF MASK"], "")

            # ESO stores CCF RV values in km/s. Convert to m/s.
            row.update({
                "object": obj,
                "target": target,
                "obs_id": obs_id,
                "date_obs": date_obs,
                "date_from_filename": date_from_name,
                "mjd_obs": mjd_obs,
                "bjd_approx": mjd_obs + 2400000.5 if np.isfinite(mjd_obs) else np.nan,
                "exptime": exptime,
                "ins_mode": ins_mode or mode_from_name,
                "mode_from_filename": mode_from_name,
                "pro_catg": pro_catg,
                "prodcatg": prodcatg,
                "rv_kms": rv,
                "rv_err_kms": rv_err,
                "berv_kms": berv,
                "ccf_fwhm_kms": fwhm,
                "ccf_fwhm_err_kms": fwhm_err,
                "ccf_bis_span_kms": bis,
                "ccf_bis_span_err_kms": bis_err,
                "rv_mps": rv * 1000 if np.isfinite(rv) else np.nan,
                "rv_err_mps": rv_err * 1000 if np.isfinite(rv_err) else np.nan,
                "berv_mps": berv * 1000 if np.isfinite(berv) else np.nan,
                "ccf_fwhm_mps": fwhm * 1000 if np.isfinite(fwhm) else np.nan,
                "ccf_fwhm_err_mps": fwhm_err * 1000 if np.isfinite(fwhm_err) else np.nan,
                "ccf_bis_span_mps": bis * 1000 if np.isfinite(bis) else np.nan,
                "ccf_bis_span_err_mps": bis_err * 1000 if np.isfinite(bis_err) else np.nan,
                "ccf_contrast": contrast,
                "ccf_contrast_err": contrast_err,
                "ccf_mask": mask,
                "n_hdu": len(hdul),
                "ok": True,
            })

    except Exception as e:
        row["error"] = repr(e)

    rows.append(row)

df = pd.DataFrame(rows).sort_values(["mjd_obs", "filename"], na_position="last")

out = OUTDIR / "espresso_ccf_header_rv.csv"
df.to_csv(out, index=False)

summary_out = METADIR / "espresso_ccf_header_rv_summary.txt"
with open(summary_out, "w") as f:
    f.write(f"Input files: {len(files)}\n")
    f.write(f"Rows: {len(df)}\n")
    f.write(f"OK: {int(df['ok'].sum())}\n")
    for c in ["rv_mps", "rv_err_mps", "berv_mps", "ccf_fwhm_mps", "ccf_contrast", "ccf_bis_span_mps", "ccf_mask"]:
        if c in df.columns:
            f.write(f"{c} finite/non-null: {df[c].notna().sum()}\n")

print("Input files:", len(files))
print("Saved:", out)
print("Saved:", summary_out)
print()
print("Finite counts:")
for c in ["rv_mps", "rv_err_mps", "berv_mps", "ccf_fwhm_mps", "ccf_contrast", "ccf_bis_span_mps", "ccf_mask"]:
    if c in df.columns:
        print(f"{c:22s}", df[c].notna().sum())

print("\nDate/MJD/RV range:")
print(df[["mjd_obs", "bjd_approx", "rv_mps", "rv_err_mps"]].agg(["min", "max", "median"]).to_string())

print("\nModes:")
for c in ["ins_mode", "pro_catg", "prodcatg", "ccf_mask"]:
    if c in df.columns:
        print("\n", c)
        print(df[c].value_counts(dropna=False).to_string())

print("\nPreview:")
cols = [
    "filename", "date_obs", "mjd_obs", "bjd_approx", "exptime",
    "ins_mode", "pro_catg", "rv_mps", "rv_err_mps", "berv_mps",
    "ccf_fwhm_mps", "ccf_contrast", "ccf_bis_span_mps", "ccf_mask"
]
cols = [c for c in cols if c in df.columns]
print(df[cols].head(20).to_string(index=False))

print("\nDONE")
