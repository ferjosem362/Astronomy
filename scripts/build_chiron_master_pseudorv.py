from pathlib import Path
import pandas as pd
import numpy as np

inp_base = Path("data_interim/qc/chiron_master_terms_qc.csv")
inp_rv   = Path("data_interim/qc/chiron_pseudorv.csv")
out      = Path("data_interim/qc/chiron_master_terms_pseudorv.csv")

base = pd.read_csv(inp_base).copy()
rvdf = pd.read_csv(inp_rv).copy()

use = rvdf[["stem","pseudo_rv_ms","pseudo_rv_err_ms"]].rename(columns={
    "pseudo_rv_ms": "rv",
    "pseudo_rv_err_ms": "rv_err",
})

df = base.drop(columns=[c for c in ["rv","rv_err"] if c in base.columns]).merge(use, on="stem", how="left")

wanted = [
    "instrument","stem","target_name","date_obs","bjd_tdb","mjd_obs",
    "rv","rv_err","berv","fwhm","contrast","bis_span","snr","airmass","drift",
    "year_sin_1","year_cos_1","year_sin_2","year_cos_2",
    "aca_sep_arcsec","aca_contam_proxy",
    "file_spectrum","file_ccf","file_bis","file_e2ds","file_drift",
    "quality_flag","notes"
]
for c in wanted:
    if c not in df.columns:
        df[c] = np.nan

df = df[wanted].sort_values("date_obs")
df.to_csv(out, index=False)

print(df[["instrument","date_obs","rv","rv_err","airmass","quality_flag"]].to_string(index=False))
print(f"\nSaved {out}")
