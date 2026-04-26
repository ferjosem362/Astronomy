from pathlib import Path
import pandas as pd
import numpy as np
import sys

YEAR = 365.25

inp = Path(sys.argv[1])
out = Path(sys.argv[2])

df = pd.read_csv(inp).copy()

for c in ["year_sin_1","year_cos_1","year_sin_2","year_cos_2"]:
    df[c] = np.nan

if "bjd_tdb" in df.columns:
    ok = df["bjd_tdb"].notna()
    if ok.any():
        t = df.loc[ok, "bjd_tdb"].astype(float).values
        x = 2.0 * np.pi * t / YEAR
        df.loc[ok, "year_sin_1"] = np.sin(x)
        df.loc[ok, "year_cos_1"] = np.cos(x)
        df.loc[ok, "year_sin_2"] = np.sin(2.0 * x)
        df.loc[ok, "year_cos_2"] = np.cos(2.0 * x)

df.to_csv(out, index=False)
print(f"saved {out}")
