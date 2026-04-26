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

fits_files = sorted(INDIR.glob("*.fits"))
print("Input CCF FITS:", len(fits_files))

HEADER_KEYS = [
    "OBJECT", "DATE-OBS", "MJD-OBS", "BJD", "EXPTIME", "RA", "DEC",
    "HIERARCH ESO OBS TARG NAME",
    "HIERARCH ESO DPR TYPE",
    "HIERARCH ESO DPR CATG",
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
    "HIERARCH ESO QC SNR AVG",
    "HIERARCH ESO OCS OBJ RV",
    "HIERARCH ESO RV START",
    "HIERARCH ESO RV STEP",
]

def hget_any(hdul, key):
    for hdu in hdul:
        try:
            if key in hdu.header:
                return hdu.header.get(key)
        except Exception:
            pass
    return np.nan

def norm_key(k):
    return re.sub(r"[^A-Za-z0-9]+", "_", k).strip("_").lower()

def extract_orig_from_filename(fn):
    # Examples:
    # ES_SCCA_3560807_2023-04-20T03:56:17.281_HR_1x1_U1.fits__ADP...
    m = re.search(r"(ES_SCCA_[^_]+_\d{4}-\d{2}-\d{2}T\d{2}[:_]\d{2}[:_]\d{2}\.\d+_[A-Z0-9]+_1x1_U\d+)", fn)
    return m.group(1).replace("_", ":", 0) if m else ""

def parse_date_from_name(fn):
    m = re.search(r"(\d{4}-\d{2}-\d{2}T\d{2}[_:]\d{2}[_:]\d{2}\.\d+)", fn)
    return m.group(1).replace("_", ":") if m else ""

rows = []
ext_rows = []
col_rows = []

for i, p in enumerate(fits_files, 1):
    row = {
        "file": str(p),
        "filename": p.name,
        "size_bytes": p.stat().st_size,
        "ok": False,
        "error": "",
        "date_from_filename": parse_date_from_name(p.name),
    }

    try:
        with fits.open(p, memmap=True) as hdul:
            row["ok"] = True
            row["n_ext"] = len(hdul)

            for key in HEADER_KEYS:
                row[norm_key(key)] = hget_any(hdul, key)

            # guardar inventario de extensiones/columnas
            for ext, hdu in enumerate(hdul):
                extname = hdu.header.get("EXTNAME", "")
                data = getattr(hdu, "data", None)
                shape = getattr(data, "shape", None)
                dtype = str(getattr(data, "dtype", "")) if data is not None else ""
                ext_rows.append({
                    "file": str(p),
                    "filename": p.name,
                    "ext": ext,
                    "extname": extname,
                    "hdu_class": type(hdu).__name__,
                    "data_shape": str(shape),
                    "data_dtype": dtype,
                })

                if hasattr(hdu, "columns"):
                    for c in hdu.columns:
                        col_rows.append({
                            "file": str(p),
                            "filename": p.name,
                            "ext": ext,
                            "extname": extname,
                            "col_name": c.name,
                            "col_format": c.format,
                            "col_unit": c.unit,
                            "col_dim": c.dim,
                        })

            # intentar extraer info de tabla si hubiera columnas CCF/RV
            table_cols = []
            for hdu in hdul:
                if hasattr(hdu, "columns") and hdu.data is not None:
                    names = list(hdu.columns.names)
                    table_cols.extend(names)
                    upper = [str(x).upper() for x in names]

                    # si hay columnas tipo RV / CCF, extraer primeras estadísticas
                    for name in names:
                        uname = str(name).upper()
                        if any(tok in uname for tok in ["RV", "CCF", "BIS", "FWHM", "CONTRAST"]):
                            try:
                                arr = np.asarray(hdu.data[name]).astype(float).ravel()
                                arr = arr[np.isfinite(arr)]
                                if arr.size:
                                    row[f"tbl_{norm_key(name)}_n"] = int(arr.size)
                                    row[f"tbl_{norm_key(name)}_median"] = float(np.nanmedian(arr))
                                    row[f"tbl_{norm_key(name)}_min"] = float(np.nanmin(arr))
                                    row[f"tbl_{norm_key(name)}_max"] = float(np.nanmax(arr))
                            except Exception:
                                pass

            row["table_columns_joined"] = "|".join(map(str, sorted(set(table_cols))))

    except Exception as e:
        row["error"] = repr(e)

    rows.append(row)

    if i % 100 == 0:
        print(f"progress {i}/{len(fits_files)}")

df = pd.DataFrame(rows)

# Convertir ESO QC RV headers a m/s si están en km/s, como en S1D.
# En tus S1D vimos RV=-24.779 km/s => -24779 m/s.
def num(s):
    return pd.to_numeric(s, errors="coerce")

mapping = {
    "rv_kms": "hierarch_eso_qc_ccf_rv",
    "rv_err_kms": "hierarch_eso_qc_ccf_rv_error",
    "berv_kms": "hierarch_eso_qc_berv",
    "ccf_fwhm_kms": "hierarch_eso_qc_ccf_fwhm",
    "ccf_fwhm_err_kms": "hierarch_eso_qc_ccf_fwhm_error",
    "ccf_bis_span_kms": "hierarch_eso_qc_ccf_bis_span",
    "ccf_bis_span_err_kms": "hierarch_eso_qc_ccf_bis_span_error",
}

for new, old in mapping.items():
    if old in df.columns:
        df[new] = num(df[old])
        df[new.replace("_kms", "_mps")] = df[new] * 1000.0

if "hierarch_eso_qc_ccf_contrast" in df.columns:
    df["ccf_contrast"] = num(df["hierarch_eso_qc_ccf_contrast"])
if "hierarch_eso_qc_ccf_contrast_error" in df.columns:
    df["ccf_contrast_err"] = num(df["hierarch_eso_qc_ccf_contrast_error"])

if "mjd_obs" in df.columns:
    df["mjd_obs"] = num(df["mjd_obs"])
    df["bjd_approx"] = df["mjd_obs"] + 2400000.5

# tabla RV limpia
rv_cols = [
    "filename", "file", "date_obs", "date_from_filename", "mjd_obs", "bjd_approx",
    "exptime", "object", "hierarch_eso_obs_targ_name", "hierarch_eso_ins_mode",
    "hierarch_eso_pro_catg", "prodcatg",
    "rv_mps", "rv_err_mps", "berv_mps",
    "ccf_fwhm_mps", "ccf_fwhm_err_mps",
    "ccf_contrast", "ccf_contrast_err",
    "ccf_bis_span_mps", "ccf_bis_span_err_mps",
    "hierarch_eso_qc_ccf_mask",
    "hierarch_eso_qc_snr", "hierarch_eso_qc_snr_avg",
    "hierarch_eso_ocs_obj_rv", "hierarch_eso_rv_start", "hierarch_eso_rv_step",
    "table_columns_joined",
]
rv_cols = [c for c in rv_cols if c in df.columns]

rv = df[df.get("rv_mps", pd.Series(index=df.index, dtype=float)).notna()].copy()
rv = rv[rv_cols].sort_values(["mjd_obs", "filename"], na_position="last")

full_out = METADIR / "espresso_ccf_fits_file_summary_full.csv"
ext_out = METADIR / "espresso_ccf_fits_extensions_full.csv"
cols_out = METADIR / "espresso_ccf_fits_columns_full.csv"
rv_out = OUTDIR / "espresso_ccf_header_rv_full.csv"

df.to_csv(full_out, index=False)
pd.DataFrame(ext_rows).to_csv(ext_out, index=False)
pd.DataFrame(col_rows).to_csv(cols_out, index=False)
rv.to_csv(rv_out, index=False)

print("\nSaved:", full_out)
print("Saved:", ext_out)
print("Saved:", cols_out)
print("Saved:", rv_out)

print("\nOK/error:")
print(df["ok"].value_counts(dropna=False).to_string())

print("\nRows with RV:", len(rv))
if len(rv):
    print("\nRV preview:")
    print(rv.head(40).to_string(index=False))

print("\nColumn inventory:")
if col_rows:
    cdf = pd.DataFrame(col_rows)
    print(cdf.groupby(["ext", "extname", "col_name"], dropna=False).size().reset_index(name="n").to_string(index=False))
else:
    print("No table columns found")

print("\nDONE")
