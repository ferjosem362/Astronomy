#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits

CCFDIR = Path("data_raw/espresso/ccf")
OUTDIR = Path("data_processed/espresso")
METADIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)
METADIR.mkdir(parents=True, exist_ok=True)

files = sorted(CCFDIR.glob("*.fits"))

HDR_KEYS = [
    "OBJECT", "DATE-OBS", "MJD-OBS", "EXPTIME", "RA", "DEC",
    "HIERARCH ESO OBS TARG NAME",
    "HIERARCH ESO INS MODE",
    "HIERARCH ESO PRO CATG",
    "PRODCATG",
    "HIERARCH ESO QC BERV",
    "HIERARCH ESO QC CCF RV",
    "HIERARCH ESO QC CCF RV ERROR",
    "HIERARCH ESO QC CCF FWHM",
    "HIERARCH ESO QC CCF FWHM ERROR",
    "HIERARCH ESO QC CCF CONTRAST",
    "HIERARCH ESO QC CCF CONTRAST ERROR",
    "HIERARCH ESO QC CCF BIS SPAN",
    "HIERARCH ESO QC CCF BIS SPAN ERROR",
    "HIERARCH ESO QC CCF MASK",
    "HIERARCH ESO QC CCF FLUX ASYMMETRY",
    "HIERARCH ESO QC CCF FLUX ASYMMETRY ERROR",
    "HIERARCH ESO QC SNR",
]

def hget(hdr, key):
    try:
        return hdr.get(key, np.nan)
    except Exception:
        return np.nan

rows = []
ext_rows = []
col_rows = []

for p in files:
    rec = {"file": str(p), "filename": p.name, "size_bytes": p.stat().st_size, "ok": False, "error": ""}
    try:
        with fits.open(p, memmap=True) as hdul:
            rec["ok"] = True
            rec["n_ext"] = len(hdul)

            # merge useful headers: primary first, then extension headers fill missing/extra
            hdr_any = {}
            for hdu in hdul:
                for k, v in hdu.header.items():
                    if k not in hdr_any:
                        hdr_any[k] = v

            for k in HDR_KEYS:
                rec[k] = hdr_any.get(k, np.nan)

            # fallback: grab any keys containing CCF/RV/BERV
            for k, v in hdr_any.items():
                ku = str(k).upper()
                if any(tok in ku for tok in ["CCF", "BERV", " RV", "RVC"]):
                    rec[f"HDR_ANY__{k}"] = v

            for i, hdu in enumerate(hdul):
                extname = hdu.header.get("EXTNAME", "")
                data = hdu.data
                er = {
                    "file": str(p),
                    "filename": p.name,
                    "ext": i,
                    "extname": extname,
                    "hdu_class": type(hdu).__name__,
                    "data_shape": str(getattr(data, "shape", "")),
                    "data_dtype": str(getattr(data, "dtype", "")),
                    "naxis": hdu.header.get("NAXIS", np.nan),
                }
                ext_rows.append(er)

                if hasattr(data, "columns"):
                    for col in data.columns:
                        col_rows.append({
                            "file": str(p),
                            "filename": p.name,
                            "ext": i,
                            "extname": extname,
                            "col_name": col.name,
                            "col_format": col.format,
                            "col_unit": col.unit,
                            "col_dim": col.dim,
                        })

    except Exception as e:
        rec["error"] = repr(e)

    rows.append(rec)

summary = pd.DataFrame(rows)
exts = pd.DataFrame(ext_rows)
cols = pd.DataFrame(col_rows)

summary_out = METADIR / "espresso_ccf_fits_file_summary.csv"
exts_out = METADIR / "espresso_ccf_fits_extensions.csv"
cols_out = METADIR / "espresso_ccf_fits_columns.csv"

summary.to_csv(summary_out, index=False)
exts.to_csv(exts_out, index=False)
cols.to_csv(cols_out, index=False)

def pick(df, *names):
    for n in names:
        if n in df.columns:
            return n
    return None

rv_col = pick(summary, "HIERARCH ESO QC CCF RV", "HDR_ANY__ESO QC CCF RV")
err_col = pick(summary, "HIERARCH ESO QC CCF RV ERROR", "HDR_ANY__ESO QC CCF RV ERROR")
berv_col = pick(summary, "HIERARCH ESO QC BERV", "HDR_ANY__ESO QC BERV")
fwhm_col = pick(summary, "HIERARCH ESO QC CCF FWHM", "HDR_ANY__ESO QC CCF FWHM")
contrast_col = pick(summary, "HIERARCH ESO QC CCF CONTRAST", "HDR_ANY__ESO QC CCF CONTRAST")
bis_col = pick(summary, "HIERARCH ESO QC CCF BIS SPAN", "HDR_ANY__ESO QC CCF BIS SPAN")
mask_col = pick(summary, "HIERARCH ESO QC CCF MASK", "HDR_ANY__ESO QC CCF MASK")

rv = pd.DataFrame()
rv["filename"] = summary["filename"]
rv["file"] = summary["file"]
rv["object"] = summary.get("OBJECT", np.nan)
rv["target"] = summary.get("HIERARCH ESO OBS TARG NAME", np.nan)
rv["date_obs"] = summary.get("DATE-OBS", np.nan)
rv["mjd_obs"] = pd.to_numeric(summary.get("MJD-OBS", np.nan), errors="coerce")
rv["exptime"] = pd.to_numeric(summary.get("EXPTIME", np.nan), errors="coerce")
rv["ins_mode"] = summary.get("HIERARCH ESO INS MODE", np.nan)
rv["pro_catg"] = summary.get("HIERARCH ESO PRO CATG", summary.get("PRODCATG", np.nan))
rv["rv_mps"] = pd.to_numeric(summary[rv_col], errors="coerce") * 1000 if rv_col else np.nan
rv["rv_err_mps"] = pd.to_numeric(summary[err_col], errors="coerce") * 1000 if err_col else np.nan
rv["berv_mps"] = pd.to_numeric(summary[berv_col], errors="coerce") * 1000 if berv_col else np.nan
rv["ccf_fwhm_mps"] = pd.to_numeric(summary[fwhm_col], errors="coerce") * 1000 if fwhm_col else np.nan
rv["ccf_contrast"] = pd.to_numeric(summary[contrast_col], errors="coerce") if contrast_col else np.nan
rv["ccf_bis_span_mps"] = pd.to_numeric(summary[bis_col], errors="coerce") * 1000 if bis_col else np.nan
rv["ccf_mask"] = summary[mask_col] if mask_col else np.nan
rv["ok"] = summary["ok"]
rv["error"] = summary["error"]

rv = rv.sort_values(["mjd_obs", "filename"], na_position="last")
rv_out = OUTDIR / "espresso_ccf_header_rv.csv"
rv.to_csv(rv_out, index=False)

print("Input CCF FITS:", len(files))
print("Saved:", summary_out)
print("Saved:", exts_out)
print("Saved:", cols_out)
print("Saved:", rv_out)

print("\nOK counts:")
print(summary["ok"].value_counts(dropna=False).to_string())

print("\nExtension names:")
if len(exts):
    print(exts["extname"].fillna("").value_counts().head(50).to_string())

print("\nColumn names:")
if len(cols):
    print(cols["col_name"].value_counts().head(80).to_string())

print("\nRV rows:")
print(rv[["filename", "date_obs", "mjd_obs", "exptime", "rv_mps", "rv_err_mps", "berv_mps", "ccf_fwhm_mps", "ccf_contrast", "ccf_bis_span_mps", "ccf_mask", "ins_mode"]].head(40).to_string(index=False))

print("\nRV non-null:", int(rv["rv_mps"].notna().sum()))
print("RV err non-null:", int(rv["rv_err_mps"].notna().sum()))
print("DONE")
