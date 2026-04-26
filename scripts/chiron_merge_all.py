from pathlib import Path
import pandas as pd

raw = pd.read_csv("data_interim/qc/chiron_master_terms_qc.csv")
bjd = pd.read_csv("data_interim/qc/chiron_master_bjd.csv")
feat = pd.read_csv("data_interim/qc/chiron_raw_features.csv")

# normaliza claves
for df in (raw, bjd, feat):
    if "stem" in df.columns:
        df["stem"] = df["stem"].astype(str)

# merge robusto
df = raw.merge(
    bjd[["stem", "mjd_obs", "bjd_tdb"]],
    on="stem",
    how="left",
    suffixes=("", "_bjd")
)

# si ya existían vacíos, rellena
if "mjd_obs_bjd" in df.columns:
    df["mjd_obs"] = df["mjd_obs"].fillna(df["mjd_obs_bjd"])
    df = df.drop(columns=["mjd_obs_bjd"])

if "bjd_tdb_bjd" in df.columns:
    df["bjd_tdb"] = df["bjd_tdb"].fillna(df["bjd_tdb_bjd"])
    df = df.drop(columns=["bjd_tdb_bjd"])

df = df.merge(feat, on=["stem", "file_spectrum"], how="left", suffixes=("", "_feat"))

out = Path("data_interim/qc/chiron_working.csv")
df = df.sort_values("date_obs").reset_index(drop=True)
df.to_csv(out, index=False)

print(df[[
    "stem","date_obs","bjd_tdb","airmass",
    "median_counts","p95_counts","max_counts","sat_frac",
    "flux_total","x_center_mass","y_center_mass","quality_flag"
]].to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
