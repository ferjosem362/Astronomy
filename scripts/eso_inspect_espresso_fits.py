#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
from astropy.io import fits
import traceback

INDIR = Path("data_raw/espresso/this")
OUTDIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)

files = sorted(INDIR.glob("*.fits"))

rows = []
ext_rows = []
col_rows = []

print("Input FITS:", len(files))

interesting_keys = [
    "OBJECT", "HIERARCH ESO OBS TARG NAME", "HIERARCH ESO DPR TYPE",
    "HIERARCH ESO DPR CATG", "HIERARCH ESO INS MODE",
    "HIERARCH ESO PRO CATG", "PRODCATG",
    "DATE-OBS", "MJD-OBS", "EXPTIME", "RA", "DEC",
    "HIERARCH ESO QC CCF RV", "HIERARCH ESO QC CCF RV ERROR",
    "HIERARCH ESO QC BERV", "HIERARCH ESO QC CCF FWHM",
    "HIERARCH ESO QC CCF CONTRAST", "HIERARCH ESO QC CCF BIS SPAN",
    "HIERARCH ESO TEL AIRM START", "HIERARCH ESO TEL AIRM END",
    "HIERARCH ESO QC SNR", "SNR",
]

for i, path in enumerate(files, 1):
    print(f"[{i:03d}/{len(files):03d}] {path.name}")

    rec = {
        "file": str(path),
        "filename": path.name,
        "size_bytes": path.stat().st_size,
        "ok": False,
        "n_ext": None,
        "error": "",
    }

    try:
        with fits.open(path, memmap=True) as hdul:
            rec["ok"] = True
            rec["n_ext"] = len(hdul)

            ph = hdul[0].header
            for k in interesting_keys:
                rec[k] = ph.get(k, None)

            # fallback: buscar claves RV/QC en cualquier extensión
            all_headers = {}
            for ext_i, hdu in enumerate(hdul):
                hdr = hdu.header
                for k in hdr.keys():
                    ku = str(k).upper()
                    if any(token in ku for token in ["RV", "BERV", "CCF", "BIS", "FWHM", "SNR"]):
                        all_headers.setdefault(k, hdr.get(k))

            for k, v in all_headers.items():
                rec[f"HDR_ANY__{k}"] = v

            for ext_i, hdu in enumerate(hdul):
                hdr = hdu.header
                data = hdu.data

                er = {
                    "file": str(path),
                    "filename": path.name,
                    "ext": ext_i,
                    "extname": hdr.get("EXTNAME", ""),
                    "hdu_class": type(hdu).__name__,
                    "data_shape": str(getattr(data, "shape", None)),
                    "data_dtype": str(getattr(data, "dtype", None)),
                    "naxis": hdr.get("NAXIS", None),
                    "pro_catg": hdr.get("HIERARCH ESO PRO CATG", hdr.get("PRODCATG", None)),
                }
                ext_rows.append(er)

                # columnas si es BinTable
                if hasattr(data, "columns"):
                    for col in data.columns:
                        col_rows.append({
                            "file": str(path),
                            "filename": path.name,
                            "ext": ext_i,
                            "extname": hdr.get("EXTNAME", ""),
                            "col_name": col.name,
                            "col_format": col.format,
                            "col_unit": getattr(col, "unit", None),
                            "col_dim": getattr(col, "dim", None),
                        })

    except Exception as e:
        rec["error"] = repr(e)
        print("  ERROR:", repr(e))
        # print(traceback.format_exc())

    rows.append(rec)

summary = pd.DataFrame(rows)
exts = pd.DataFrame(ext_rows)
cols = pd.DataFrame(col_rows)

summary_out = OUTDIR / "espresso_fits_file_summary.csv"
exts_out = OUTDIR / "espresso_fits_extensions.csv"
cols_out = OUTDIR / "espresso_fits_columns.csv"

summary.to_csv(summary_out, index=False)
exts.to_csv(exts_out, index=False)
cols.to_csv(cols_out, index=False)

print("\nSaved:", summary_out)
print("Saved:", exts_out)
print("Saved:", cols_out)

print("\nOK counts:")
print(summary["ok"].value_counts(dropna=False).to_string())

print("\nExtension names:")
if len(exts):
    print(exts.groupby(["ext", "extname", "hdu_class"], dropna=False).size().reset_index(name="n").to_string(index=False))

print("\nColumn names:")
if len(cols):
    print(cols.groupby(["ext", "extname", "col_name"], dropna=False).size().reset_index(name="n").sort_values(["ext", "extname", "col_name"]).to_string(index=False))

print("\nLikely RV/header columns present:")
for c in summary.columns:
    cu = c.upper()
    if any(t in cu for t in ["RV", "BERV", "CCF", "BIS", "FWHM", "SNR"]):
        nn = summary[c].notna().sum()
        if nn:
            print(c, "non-null:", nn, "examples:", summary[c].dropna().head(5).tolist())

print("\nDONE")
