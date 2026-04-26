from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/harps_aux_master_with_bis.csv")
out = Path("data_interim/qc/harps_master.csv")

df = pd.read_csv(inp)

df["instrument"] = "HARPS"
df["target_name"] = df["object"]
df["bjd_tdb"] = df["bjd"]
df["drift"] = np.nan
df["aca_sep_arcsec"] = np.nan
df["aca_contam_proxy"] = np.nan
df["quality_flag"] = 0
df["notes"] = ""

keep = [
    "instrument","stem","target_name","date_obs","bjd_tdb","mjd_obs",
    "rv","rv_err","berv","fwhm","contrast","bis_span","snr","airmass","drift",
    "aca_sep_arcsec","aca_contam_proxy",
    "ccf_file","bis_file","s1d_file","e2ds_file","quality_flag","notes"
]

if "bis_span" not in df.columns:
    if "bis_rv" in df.columns:
        df["bis_span"] = np.nan
    else:
        df["bis_span"] = np.nan

rename = {
    "ccf_file":"file_ccf",
    "bis_file":"file_bis",
    "s1d_file":"file_spectrum",
    "e2ds_file":"file_e2ds",
}
df = df.rename(columns=rename)

keep = [
    "instrument","stem","target_name","date_obs","bjd_tdb","mjd_obs",
    "rv","rv_err","berv","fwhm","contrast","bis_span","snr","airmass","drift",
    "aca_sep_arcsec","aca_contam_proxy",
    "file_spectrum","file_ccf","file_bis","file_e2ds","quality_flag","notes"
]

df[keep].sort_values("bjd_tdb").to_csv(out, index=False)
print(f"saved {out}")
