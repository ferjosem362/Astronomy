from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_raw/metadata/chiron_rvs.csv")
out = Path("data_interim/qc/chiron_master.csv")

df = pd.read_csv(inp)

# espera columnas razonables: bjd_tdb, rv, rv_err, snr opcionales
if "bjd_tdb" not in df.columns:
    raise SystemExit("chiron_rvs.csv needs at least bjd_tdb, rv, rv_err")

outdf = pd.DataFrame()
outdf["instrument"] = "CHIRON"
outdf["stem"] = df["bjd_tdb"].astype(str)
outdf["target_name"] = df["target_name"] if "target_name" in df.columns else "ALPCENB"
outdf["date_obs"] = df["date_obs"] if "date_obs" in df.columns else np.nan
outdf["bjd_tdb"] = df["bjd_tdb"]
outdf["mjd_obs"] = df["mjd_obs"] if "mjd_obs" in df.columns else np.nan
outdf["rv"] = df["rv"]
outdf["rv_err"] = df["rv_err"]
outdf["berv"] = df["berv"] if "berv" in df.columns else np.nan
outdf["fwhm"] = df["fwhm"] if "fwhm" in df.columns else np.nan
outdf["contrast"] = df["contrast"] if "contrast" in df.columns else np.nan
outdf["bis_span"] = df["bis_span"] if "bis_span" in df.columns else np.nan
outdf["snr"] = df["snr"] if "snr" in df.columns else np.nan
outdf["airmass"] = df["airmass"] if "airmass" in df.columns else np.nan
outdf["drift"] = np.nan
outdf["aca_sep_arcsec"] = np.nan
outdf["aca_contam_proxy"] = np.nan
outdf["file_spectrum"] = df["file_spectrum"] if "file_spectrum" in df.columns else np.nan
outdf["file_ccf"] = np.nan
outdf["file_bis"] = np.nan
outdf["file_e2ds"] = np.nan
outdf["file_drift"] = np.nan
outdf["quality_flag"] = 0
outdf["notes"] = ""

out.parent.mkdir(parents=True, exist_ok=True)
outdf.to_csv(out, index=False)
print(f"saved {out} rows={len(outdf)}")
print(outdf.head(20).to_string(index=False))
