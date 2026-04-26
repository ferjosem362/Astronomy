#!/usr/bin/env python3
from pathlib import Path
import re
import numpy as np
import pandas as pd
from astropy.io import fits

IN_GLOB = "data_raw/espresso/this/*.fits"
OUTDIR = Path("data_processed/espresso")
METADIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)
METADIR.mkdir(parents=True, exist_ok=True)

def hget(h, key, default=np.nan):
    return h[key] if key in h else default

def product_hint(filename, hdr):
    fn = Path(filename).name
    m = re.search(r"_(ES_[A-Z0-9]+)_", fn)
    if m:
        return m.group(1)
    pc = str(hget(hdr, "HIERARCH ESO PRO CATG", ""))
    if pc:
        return pc
    return "UNKNOWN"

def as_float(x):
    try:
        if x is None:
            return np.nan
        return float(x)
    except Exception:
        return np.nan

rows = []
rv_rows = []

files = sorted(Path("data_raw/espresso/this").glob("*.fits"))
print("Input FITS:", len(files))

for fp in files:
    try:
        with fits.open(fp, memmap=True) as hdul:
            hdr = hdul[0].header
            data = hdul[1].data if len(hdul) > 1 else None

            hint = product_hint(fp.name, hdr)
            row = {
                "file": str(fp),
                "filename": fp.name,
                "product_hint": hint,
                "pro_catg": hget(hdr, "HIERARCH ESO PRO CATG", ""),
                "object": hget(hdr, "OBJECT", ""),
                "date_obs": hget(hdr, "DATE-OBS", ""),
                "mjd_obs": as_float(hget(hdr, "MJD-OBS", np.nan)),
                "exptime": as_float(hget(hdr, "EXPTIME", np.nan)),
                "ra": as_float(hget(hdr, "RA", np.nan)),
                "dec": as_float(hget(hdr, "DEC", np.nan)),
                "ins_mode": hget(hdr, "HIERARCH ESO INS MODE", ""),
                "snr_header": as_float(hget(hdr, "SNR", hget(hdr, "HIERARCH ESO QC SNR", np.nan))),
                "berv_kms": as_float(hget(hdr, "HIERARCH ESO QC BERV", np.nan)),
                "ccf_rv_kms": as_float(hget(hdr, "HIERARCH ESO QC CCF RV", np.nan)),
                "ccf_rv_err_kms": as_float(hget(hdr, "HIERARCH ESO QC CCF RV ERROR", np.nan)),
                "ccf_fwhm_kms": as_float(hget(hdr, "HIERARCH ESO QC CCF FWHM", np.nan)),
                "ccf_contrast": as_float(hget(hdr, "HIERARCH ESO QC CCF CONTRAST", np.nan)),
                "ccf_bis_span_kms": as_float(hget(hdr, "HIERARCH ESO QC CCF BIS SPAN", np.nan)),
            }

            if data is not None:
                names = list(data.names or [])
                row["columns"] = ",".join(names)

                wave = np.asarray(data["WAVE"][0], dtype=float) if "WAVE" in names else np.array([])
                flux = np.asarray(data["FLUX"][0], dtype=float) if "FLUX" in names else np.array([])
                err  = np.asarray(data["ERR"][0], dtype=float)  if "ERR"  in names else np.array([])
                snr  = np.asarray(data["SNR"][0], dtype=float)  if "SNR"  in names else np.array([])
                qual = np.asarray(data["QUAL"][0])              if "QUAL" in names else np.array([])

                finite_wave = np.isfinite(wave)
                finite_flux = np.isfinite(flux)
                good = finite_wave & finite_flux
                if qual.size == wave.size:
                    good = good & (qual == 0)

                row.update({
                    "n_pix": int(wave.size),
                    "wave_min": float(np.nanmin(wave)) if wave.size else np.nan,
                    "wave_max": float(np.nanmax(wave)) if wave.size else np.nan,
                    "flux_finite_frac": float(finite_flux.mean()) if flux.size else np.nan,
                    "qual0_frac": float((qual == 0).mean()) if qual.size else np.nan,
                    "snr_median_array": float(np.nanmedian(snr)) if snr.size else np.nan,
                    "flux_median_good": float(np.nanmedian(flux[good])) if good.any() else np.nan,
                    "err_median_good": float(np.nanmedian(err[good])) if err.size == wave.size and good.any() else np.nan,
                })

            rows.append(row)

            if np.isfinite(row["ccf_rv_kms"]):
                rv = row.copy()
                rv["rv_mps"] = row["ccf_rv_kms"] * 1000.0
                rv["rv_err_mps"] = row["ccf_rv_err_kms"] * 1000.0
                rv["berv_mps"] = row["berv_kms"] * 1000.0 if np.isfinite(row["berv_kms"]) else np.nan
                rv["ccf_fwhm_mps"] = row["ccf_fwhm_kms"] * 1000.0 if np.isfinite(row["ccf_fwhm_kms"]) else np.nan
                rv["ccf_bis_span_mps"] = row["ccf_bis_span_kms"] * 1000.0 if np.isfinite(row["ccf_bis_span_kms"]) else np.nan
                rv_rows.append(rv)

    except Exception as e:
        rows.append({
            "file": str(fp),
            "filename": fp.name,
            "error": repr(e),
        })

index = pd.DataFrame(rows).sort_values(["mjd_obs", "filename"], na_position="last")
rvtab = pd.DataFrame(rv_rows).sort_values(["mjd_obs", "filename"], na_position="last")

index_out = METADIR / "espresso_s1d_index.csv"
rv_out = OUTDIR / "espresso_pipeline_header_rv.csv"

index.to_csv(index_out, index=False)
rvtab.to_csv(rv_out, index=False)

print("Saved:", index_out, len(index))
print("Saved:", rv_out, len(rvtab))

print("\nProduct/pro_catg counts:")
print(index.groupby(["product_hint", "pro_catg"], dropna=False).size().reset_index(name="n").to_string(index=False))

print("\nS1D index compact:")
cols = [
    "filename", "product_hint", "pro_catg", "date_obs", "mjd_obs", "exptime",
    "ins_mode", "n_pix", "wave_min", "wave_max", "snr_header",
    "snr_median_array", "qual0_frac", "flux_finite_frac"
]
cols = [c for c in cols if c in index.columns]
print(index[cols].head(30).to_string(index=False))

print("\nPipeline RV rows:")
rvcols = [
    "filename", "product_hint", "pro_catg", "date_obs", "mjd_obs", "exptime",
    "rv_mps", "rv_err_mps", "berv_mps", "ccf_fwhm_mps",
    "ccf_contrast", "ccf_bis_span_mps", "ins_mode", "snr_header"
]
rvcols = [c for c in rvcols if c in rvtab.columns]
if len(rvtab):
    print(rvtab[rvcols].to_string(index=False))
else:
    print("No RV rows found.")

print("\nDONE")
