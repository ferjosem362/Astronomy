#!/usr/bin/env python3
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import urllib.request
import pandas as pd

MAN = Path("data_raw/metadata/eso_espresso_ccf_manifest.csv")
OUT = Path("data_raw/metadata/eso_espresso_ccf_download_status.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(MAN)

def download_one(row):
    dest = Path(row["local_path"])
    dest.parent.mkdir(parents=True, exist_ok=True)
    expected = row.get("datalink_content_length", None)
    try:
        expected = int(expected) if pd.notna(expected) else None
    except Exception:
        expected = None

    if dest.exists():
        sz = dest.stat().st_size
        if expected is None or sz == expected:
            return {**row.to_dict(), "download_status": "exists_ok", "bytes_got": sz, "download_error": ""}

    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err = ""
    for attempt in range(1, 5):
        try:
            req = urllib.request.Request(
                row["download_url"],
                headers={"User-Agent": "alpha-cen-b-reanalysis/0.1"}
            )
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1024 * 512)
                    if not chunk:
                        break
                    f.write(chunk)

            tmp.replace(dest)
            sz = dest.stat().st_size
            if expected is not None and sz != expected:
                last_err = f"size_mismatch got={sz} expected={expected}"
                time.sleep(attempt)
                continue

            return {**row.to_dict(), "download_status": "downloaded_ok", "bytes_got": sz, "download_error": ""}
        except Exception as e:
            last_err = repr(e)
            time.sleep(min(2 * attempt, 10))

    sz = dest.stat().st_size if dest.exists() else 0
    return {**row.to_dict(), "download_status": "ERROR", "bytes_got": sz, "download_error": last_err}

rows = []
workers = 8
print("Input rows:", len(df))
print("Workers:", workers)

with ThreadPoolExecutor(max_workers=workers) as ex:
    futs = [ex.submit(download_one, row) for _, row in df.iterrows()]
    for n, fut in enumerate(as_completed(futs), 1):
        res = fut.result()
        rows.append(res)
        if n % 50 == 0 or res["download_status"] == "ERROR":
            print(f"progress {n}/{len(df)} {res['download_status']}")

status = pd.DataFrame(rows)
status.to_csv(OUT, index=False)

print("Saved:", OUT)
print(status["download_status"].value_counts(dropna=False).to_string())
print("Downloaded GB:", round(status["bytes_got"].fillna(0).sum() / 1024**3, 3))
bad = status[~status["download_status"].isin(["downloaded_ok", "exists_ok"])]
print("Bad rows:", len(bad))
if len(bad):
    print(bad[["download_url", "local_path", "download_status", "download_error"]].head(20).to_string(index=False))
print("DONE")
