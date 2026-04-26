#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import pyvo

OUTDIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)

TAP_URL = "https://archive.eso.org/tap_obs"
svc = pyvo.dal.TAPService(TAP_URL)

RA = 219.902066
DEC = -60.833976
RADIUS_DEG = 0.08   # un poco más ancho para no perder ESPRESSO por coordenadas/header

QUERY = f"""
SELECT TOP 10000
    obs_publisher_did,
    obs_id,
    obs_title,
    dataproduct_type,
    calib_level,
    target_name,
    s_ra,
    s_dec,
    s_fov,
    em_min,
    em_max,
    em_res_power,
    t_min,
    t_max,
    t_exptime,
    instrument_name,
    facility_name,
    obs_collection,
    o_ucd,
    access_url,
    access_format,
    access_estsize
FROM ivoa.ObsCore
WHERE
    1 = CONTAINS(
        POINT('ICRS', s_ra, s_dec),
        CIRCLE('ICRS', {RA}, {DEC}, {RADIUS_DEG})
    )
    AND dataproduct_type = 'spectrum'
    AND (
        instrument_name = 'ESPRESSO'
        OR obs_collection = 'ESPRESSO'
        OR instrument_name LIKE '%ESPRESSO%'
        OR obs_collection LIKE '%ESPRESSO%'
        OR obs_title LIKE '%ES_SOB%'
        OR obs_title LIKE '%ES_SFLX%'
    )
ORDER BY t_min
"""

print("Querying ESO TAP for ESPRESSO around Alpha Cen B...")
print(QUERY)

res = svc.search(QUERY)
df = res.to_table().to_pandas()

out_csv = OUTDIR / "eso_espresso_hd128621_obscore.csv"
df.to_csv(out_csv, index=False)

print("\nRows:", len(df))
print("Saved:", out_csv)

pd.set_option("display.max_colwidth", 180)
pd.set_option("display.max_rows", 200)

for c in [
    "obs_collection", "instrument_name", "target_name",
    "dataproduct_type", "calib_level", "em_res_power",
    "t_min", "t_exptime", "access_estsize", "obs_title", "access_url"
]:
    if c not in df.columns:
        df[c] = None

print("\nCounts by collection/instrument/type/level:")
if len(df):
    print(
        df.groupby(
            ["obs_collection", "instrument_name", "dataproduct_type", "calib_level"],
            dropna=False
        )
        .size()
        .reset_index(name="n")
        .sort_values("n", ascending=False)
        .to_string(index=False)
    )

print("\nCounts target_name:")
if len(df):
    print(df["target_name"].value_counts(dropna=False).head(50).to_string())

cols = [
    "obs_collection", "instrument_name", "target_name",
    "dataproduct_type", "calib_level", "em_res_power",
    "t_min", "t_exptime", "access_estsize", "obs_title", "access_url"
]

print("\nPreview:")
if len(df):
    print(df[cols].head(120).to_string(index=False))

print("\nRows with likely S2D/S2DA/SOBF/SFLX hints:")
if len(df):
    mask = (
        df["obs_title"].fillna("").str.contains("S2D|S2DA|SOBF|SFLX", case=False, regex=True)
        | df["obs_id"].fillna("").str.contains("S2D|S2DA|SOBF|SFLX", case=False, regex=True)
        | df["access_url"].fillna("").str.contains("S2D|S2DA|SOBF|SFLX", case=False, regex=True)
    )
    print(df.loc[mask, cols].head(200).to_string(index=False))

print("\nDONE")
