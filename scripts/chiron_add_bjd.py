from pathlib import Path
import pandas as pd
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

inp = Path("data_interim/qc/chiron_master_raw.csv")
out = Path("data_interim/qc/chiron_master_bjd.csv")

df = pd.read_csv(inp).copy()

# CTIO 1.5m / Cerro Tololo aproximado
site = EarthLocation.from_geodetic(lon=-70.815*u.deg, lat=-30.165*u.deg, height=2215*u.m)

# Alpha Cen B aproximado
target = SkyCoord(ra="14h39m35.08s", dec="-60d50m13.8s", frame="icrs")

mjd_vals = []
bjd_vals = []

for _, row in df.iterrows():
    date_obs = row.get("date_obs", None)
    if pd.isna(date_obs):
        mjd_vals.append(np.nan)
        bjd_vals.append(np.nan)
        continue
    try:
        t = Time(str(date_obs), format="isot", scale="utc", location=site)
        mjd_vals.append(float(t.utc.mjd))
        ltt = t.light_travel_time(target)
        bjd_tdb = t.tdb + ltt
        bjd_vals.append(float(bjd_tdb.jd))
    except Exception:
        mjd_vals.append(np.nan)
        bjd_vals.append(np.nan)

df["mjd_obs"] = mjd_vals
df["bjd_tdb"] = bjd_vals

df.to_csv(out, index=False)
print(df[["instrument","target_name","date_obs","mjd_obs","bjd_tdb"]].head(10).to_string(index=False))
print(f"\nSaved {out}")
