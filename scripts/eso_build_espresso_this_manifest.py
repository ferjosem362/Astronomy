#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import re

INCSV = Path("data_raw/metadata/eso_espresso_hd128621_datalinks_expanded.csv")
OUTDIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INCSV)

m = (
    df["datalink_semantics"].astype(str).eq("#this")
    & df["datalink_description"].astype(str).str.contains("Requested file", case=False, na=False)
    & df["datalink_content_type"].astype(str).str.contains("fits", case=False, na=False)
    & df["datalink_access_url"].astype(str).str.startswith(("http://", "https://"))
)

man = df[m].drop_duplicates(subset=["datalink_access_url"]).copy()

def extract_adp(url):
    mm = re.search(r"(ADP\.[A-Za-z0-9T:\.\-]+)", str(url))
    return mm.group(1) if mm else ""

def safe(s):
    s = str(s).strip()
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s.strip("_")

man["adp_id"] = man["datalink_access_url"].map(extract_adp)
man["espresso_product_hint"] = "UNKNOWN"

title = man["obs_title"].fillna("").astype(str)
man.loc[title.str.contains("ES_SOBF", case=False, na=False), "espresso_product_hint"] = "ES_SOBF"
man.loc[title.str.contains("ES_SFLX", case=False, na=False), "espresso_product_hint"] = "ES_SFLX"

def local_filename(row):
    title = safe(row.get("obs_title", ""))
    adp = safe(row.get("adp_id", "")).replace(":", "-")
    if title and title.lower() != "nan":
        return f"{title}__{adp}.fits"
    return f"{adp}.fits"

man["download_url"] = man["datalink_access_url"]
man["local_filename"] = man.apply(local_filename, axis=1)
man["local_path"] = "data_raw/espresso/this/" + man["local_filename"]

keep = [
    "source_index",
    "target_name",
    "obs_id",
    "obs_title",
    "espresso_product_hint",
    "t_min",
    "t_max",
    "t_exptime",
    "em_res_power",
    "adp_id",
    "download_url",
    "local_filename",
    "local_path",
    "datalink_content_type",
    "datalink_content_length",
]
keep = [c for c in keep if c in man.columns]

man = man[keep].sort_values(["t_min", "obs_title"], na_position="last").copy()

out = OUTDIR / "eso_espresso_hd128621_this_manifest.csv"
man.to_csv(out, index=False)

print("Input rows:", len(df))
print("Manifest #this rows:", len(man))
print("Saved:", out)

print("\nCounts product hint:")
print(man["espresso_product_hint"].value_counts(dropna=False).to_string())

print("\nCounts target:")
print(man["target_name"].value_counts(dropna=False).to_string())

print("\nTotal size GB approx:")
print(round(man["datalink_content_length"].fillna(0).sum() / 1024**3, 3))

pd.set_option("display.max_colwidth", 180)
print("\nPreview:")
print(man.head(40).to_string(index=False))

print("\nDONE")
