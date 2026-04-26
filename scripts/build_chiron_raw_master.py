from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/chiron_inventory.csv")
out = Path("data_interim/qc/chiron_master_raw.csv")

df = pd.read_csv(inp)
outdf = pd.DataFrame(index=df.index)

outdf["instrument"] = "CHIRON"
outdf["stem"] = df["filename"].astype(str).str.replace(".fits.fz", "", regex=False)
outdf["target_name"] = df["object0"].fillna(df["object1"])
outdf["date_obs"] = df["date_obs0"].fillna(df["date_obs1"])
outdf["bjd_tdb"] = np.nan
outdf["mjd_obs"] = df["mjd_obs0"].fillna(df["mjd_obs1"])
outdf["rv"] = np.nan
outdf["rv_err"] = np.nan
outdf["berv"] = np.nan
outdf["fwhm"] = np.nan
outdf["contrast"] = np.nan
outdf["bis_span"] = np.nan
outdf["snr"] = np.nan
outdf["airmass"] = df["airmass0"].fillna(df["airmass1"])
outdf["drift"] = np.nan
outdf["aca_sep_arcsec"] = np.nan
outdf["aca_contam_proxy"] = np.nan
outdf["file_spectrum"] = df["file"]
outdf["file_ccf"] = np.nan
outdf["file_bis"] = np.nan
outdf["file_e2ds"] = np.nan
outdf["file_drift"] = np.nan
outdf["quality_flag"] = 0
outdf["notes"] = "CHIRON raw only; RV not yet extracted"

out.parent.mkdir(parents=True, exist_ok=True)
outdf.to_csv(out, index=False)
print(outdf.head(10).to_string(index=False))
print(f"\nSaved {out} rows={len(outdf)}")
