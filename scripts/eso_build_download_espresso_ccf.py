#!/usr/bin/env python3
from pathlib import Path
import hashlib
import re
import time
import urllib.request
import urllib.error

import pandas as pd

INV = Path("data_raw/metadata/eso_espresso_auxiliary_inventory.csv")
OUTDIR = Path("data_raw/espresso/ccf")
METADIR = Path("data_raw/metadata")
LOGDIR = Path("logs")

OUTDIR.mkdir(parents=True, exist_ok=True)
METADIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INV)

m = (
    df["datalink_eso_category"].astype(str).eq("ANCILLARY.CCF")
    & df["datalink_content_type"].astype(str).str.contains("fits", case=False, na=False)
    & df["datalink_access_url"].astype(str).str.startswith(("http://", "https://"))
)

ccf = df[m].drop_duplicates(subset=["datalink_access_url"]).copy()

def safe(s):
    s = str(s).strip()
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", s)
    return s.strip("_")

def extract_file_id(url):
    mm = re.search(r"(ADP\.[A-Za-z0-9T:\.\-]+)", str(url))
    return mm.group(1) if mm else hashlib.sha1(str(url).encode()).hexdigest()[:16]

ccf["file_id"] = ccf["datalink_access_url"].map(extract_file_id)
ccf["download_url"] = ccf["datalink_access_url"]

def local_filename(row):
    orig = safe(row.get("datalink_eso_origfile", ""))
    fid = safe(row.get("file_id", "")).replace(":", "-")
    if orig and orig.lower() != "nan":
        return f"{orig.replace('.fits','')}__{fid}.fits"
    return f"{fid}.fits"

ccf["local_filename"] = ccf.apply(local_filename, axis=1)
ccf["local_path"] = ccf["local_filename"].map(lambda x: str(OUTDIR / x))

keep = [
    "source_index", "target_name", "obs_id", "obs_title",
    "t_min", "t_max", "t_exptime", "em_res_power",
    "datalink_eso_origfile", "datalink_eso_category",
    "datalink_content_type", "datalink_content_length",
    "file_id", "download_url", "local_filename", "local_path",
]
keep = [c for c in keep if c in ccf.columns]
ccf = ccf[keep].sort_values(["t_min", "datalink_eso_origfile"], na_position="last").copy()

manifest_path = METADIR / "eso_espresso_ccf_manifest.csv"
ccf.to_csv(manifest_path, index=False)

print("CCF manifest rows:", len(ccf))
print("Saved:", manifest_path)
print("Approx GB:", round(ccf["datalink_content_length"].fillna(0).sum() / 1024**3, 3))
print("First rows:")
pd.set_option("display.max_colwidth", 160)
print(ccf.head(20).to_string(index=False))

def download_one(url, dest, expected_size=None, timeout=120, retries=4):
    dest = Path(dest)
    tmp = dest.with_suffix(dest.suffix + ".part")

    if dest.exists():
        if expected_size and dest.stat().st_size == int(expected_size):
            return "exists_ok", dest.stat().st_size, ""
        if dest.stat().st_size > 0 and not expected_size:
            return "exists_present", dest.stat().st_size, ""

    last_err = ""
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "alpha-cen-b-reanalysis/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1024 * 512)
                    if not chunk:
                        break
                    f.write(chunk)
            tmp.replace(dest)

            sz = dest.stat().st_size
            if expected_size and sz != int(expected_size):
                last_err = f"size_mismatch got={sz} expected={int(expected_size)}"
                continue
            return "downloaded_ok", sz, ""
        except Exception as e:
            last_err = repr(e)
            time.sleep(min(2 * attempt, 10))

    return "ERROR", dest.stat().st_size if dest.exists() else 0, last_err

rows = []
for i, row in ccf.iterrows():
    status, bytes_got, err = download_one(
        row["download_url"],
        row["local_path"],
        row.get("datalink_content_length", None),
    )
    rows.append({
        **row.to_dict(),
        "download_status": status,
        "bytes_got": bytes_got,
        "download_error": err,
    })
    if len(rows) % 100 == 0:
        print(f"progress {len(rows)}/{len(ccf)} latest={status}")

status_df = pd.DataFrame(rows)
status_path = METADIR / "eso_espresso_ccf_download_status.csv"
status_df.to_csv(status_path, index=False)

print("\nSaved:", status_path)
print(status_df["download_status"].value_counts(dropna=False).to_string())
print("Downloaded GB:", round(status_df["bytes_got"].fillna(0).sum() / 1024**3, 3))
bad = status_df[~status_df["download_status"].isin(["downloaded_ok", "exists_ok", "exists_present"])]
print("Bad rows:", len(bad))
if len(bad):
    print(bad[["download_url", "local_path", "download_status", "download_error"]].head(20).to_string(index=False))

print("DONE")
