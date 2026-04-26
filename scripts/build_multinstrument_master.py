from pathlib import Path
import pandas as pd

files = [
    Path("data_interim/qc/harps_master_terms.csv"),
    Path("data_interim/qc/espresso_master_terms.csv"),
    Path("data_interim/qc/uves_master_terms.csv"),
    Path("data_interim/qc/chiron_master_terms.csv"),
]

dfs = []
for fp in files:
    if fp.exists():
        dfs.append(pd.read_csv(fp))
        print("using", fp)

if not dfs:
    raise SystemExit("no master terms files found")

df = pd.concat(dfs, ignore_index=True, sort=False)

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
        df[c] = pd.NA

df = df[wanted].sort_values(["instrument","bjd_tdb","date_obs"], na_position="last")

out = Path("data_interim/qc/multinstrument_master.csv")
df.to_csv(out, index=False)

print(f"saved {out} rows={len(df)}")
print(df.groupby("instrument").size())
