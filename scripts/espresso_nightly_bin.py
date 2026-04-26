from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/multinstrument_rvonly_centered.csv")
out = Path("data_interim/qc/espresso_nightly_binned.csv")

df = pd.read_csv(inp).copy()
df = df[df["instrument"] == "ESPRESSO"].copy()
df = df[df["bjd_tdb"].notna() & df["rv_centered"].notna() & df["rv_err"].notna()].copy()

df["night"] = np.floor(df["mjd_obs"]).astype(int)

rows = []
for night, g in df.groupby("night"):
    w = 1.0 / (g["rv_err"].values ** 2)
    rvw = np.sum(w * g["rv_centered"].values) / np.sum(w)
    errw = (1.0 / np.sum(w)) ** 0.5
    bjdw = np.sum(w * g["bjd_tdb"].values) / np.sum(w)

    rows.append({
        "instrument": "ESPRESSO",
        "night": night,
        "n_exp": len(g),
        "date_obs_min": g["date_obs"].min(),
        "date_obs_max": g["date_obs"].max(),
        "bjd_tdb": bjdw,
        "rv": rvw,
        "rv_err": errw,
        "snr_med": g["snr"].median()
    })

outdf = pd.DataFrame(rows).sort_values("bjd_tdb")
outdf.to_csv(out, index=False)

print(outdf.to_string(index=False))
print(f"\nSaved {out} rows={len(outdf)}")
