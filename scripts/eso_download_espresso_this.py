#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import requests
import time
import hashlib
import os
import sys

MANIFEST = Path("data_raw/metadata/eso_espresso_hd128621_this_manifest.csv")
OUTDIR = Path("data_raw/espresso/this")
LOGDIR = Path("logs")
OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

STATUS_CSV = Path("data_raw/metadata/eso_espresso_hd128621_this_download_status.csv")

TIMEOUT = 120
CHUNK = 1024 * 1024
SLEEP = 0.5
MAX_RETRIES = 4

df = pd.read_csv(MANIFEST)

required = ["download_url", "local_path", "local_filename"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise RuntimeError(f"Missing columns in manifest: {missing}")

session = requests.Session()
session.headers.update({
    "User-Agent": "alpha_cen_b_reanalysis/0.1"
})

def file_size(path):
    try:
        return path.stat().st_size
    except FileNotFoundError:
        return 0

def download_one(url, outpath, expected_size=None):
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    existing = file_size(outpath)

    if expected_size and existing == int(expected_size):
        return "exists_ok", existing, ""

    # Si hay parcial, intentamos Range resume
    headers = {}
    mode = "wb"
    if existing > 0:
        headers["Range"] = f"bytes={existing}-"
        mode = "ab"

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with session.get(url, headers=headers, stream=True, timeout=TIMEOUT, allow_redirects=True) as r:
                # Si pedimos Range y el servidor no lo acepta, recomenzamos limpio
                if existing > 0 and r.status_code == 200:
                    mode = "wb"
                    existing = 0

                if r.status_code not in (200, 206):
                    return "http_error", file_size(outpath), f"HTTP {r.status_code}"

                with open(outpath, mode) as f:
                    for chunk in r.iter_content(chunk_size=CHUNK):
                        if chunk:
                            f.write(chunk)

            final_size = file_size(outpath)

            if expected_size and final_size != int(expected_size):
                # Quizá hubo corte. Reintentamos con Range desde el tamaño actual.
                existing = final_size
                headers = {"Range": f"bytes={existing}-"}
                mode = "ab"
                time.sleep(SLEEP * attempt)
                continue

            return "downloaded_ok", final_size, ""

        except Exception as e:
            err = repr(e)
            existing = file_size(outpath)
            headers = {"Range": f"bytes={existing}-"} if existing > 0 else {}
            mode = "ab" if existing > 0 else "wb"
            time.sleep(SLEEP * attempt)

    return "failed", file_size(outpath), err if "err" in locals() else "unknown error"

rows = []
n = len(df)

print(f"Input manifest rows: {n}")
print(f"Output dir: {OUTDIR}")

for i, row in df.iterrows():
    url = str(row["download_url"])
    outpath = Path(str(row["local_path"]))
    expected = row.get("datalink_content_length", None)

    if pd.isna(expected):
        expected = None
    else:
        expected = int(float(expected))

    label = row.get("obs_title", "")
    if pd.isna(label) or not str(label).strip():
        label = row.get("adp_id", outpath.name)

    print(f"[{i+1:03d}/{n:03d}] {label}")
    print(f"  -> {outpath}")

    status, size, err = download_one(url, outpath, expected)

    print(f"  status={status} size={size} expected={expected} err={err}")

    rec = row.to_dict()
    rec.update({
        "download_status": status,
        "downloaded_size": size,
        "expected_size": expected,
        "download_error": err,
        "exists": outpath.exists(),
    })
    rows.append(rec)

    pd.DataFrame(rows).to_csv(STATUS_CSV, index=False)
    time.sleep(SLEEP)

status = pd.DataFrame(rows)
status.to_csv(STATUS_CSV, index=False)

print("\nSaved:", STATUS_CSV)
print("\nStatus counts:")
print(status["download_status"].value_counts(dropna=False).to_string())

print("\nDownloaded bytes:")
print(int(status["downloaded_size"].fillna(0).sum()))

print("\nDownloaded GB:")
print(round(status["downloaded_size"].fillna(0).sum() / 1024**3, 3))

bad = status[~status["download_status"].isin(["exists_ok", "downloaded_ok"])]
print("\nBad rows:", len(bad))
if len(bad):
    cols = [c for c in ["obs_title", "adp_id", "download_status", "downloaded_size", "expected_size", "download_error", "download_url"] if c in bad.columns]
    print(bad[cols].to_string(index=False))

print("\nDONE")
