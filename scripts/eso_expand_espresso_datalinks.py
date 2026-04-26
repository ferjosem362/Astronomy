#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import requests
from astropy.io.votable import parse_single_table
from io import BytesIO

INCSV = Path("data_raw/metadata/eso_espresso_hd128621_obscore.csv")
OUTDIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)

if not INCSV.exists():
    raise FileNotFoundError(f"Missing {INCSV}. Run eso_query_espresso_products.py first.")

df = pd.read_csv(INCSV)

rows = []
print(f"Input rows: {len(df)}")

for i, r in df.iterrows():
    access_url = str(r.get("access_url", "")).strip()
    obs_title = str(r.get("obs_title", "")).strip()
    obs_id = str(r.get("obs_id", "")).strip()
    target_name = str(r.get("target_name", "")).strip()

    if not access_url or access_url.lower() == "nan":
        continue

    print(f"[{i+1}/{len(df)}] {obs_title or obs_id}")

    try:
        resp = requests.get(access_url, timeout=60)
        resp.raise_for_status()

        table = parse_single_table(BytesIO(resp.content)).to_table()
        dl = table.to_pandas()

        for _, d in dl.iterrows():
            row = {
                "source_index": i,
                "target_name": target_name,
                "obs_id": obs_id,
                "obs_title": obs_title,
                "t_min": r.get("t_min"),
                "t_max": r.get("t_max"),
                "t_exptime": r.get("t_exptime"),
                "em_res_power": r.get("em_res_power"),
                "obscore_access_url": access_url,
            }

            for c in dl.columns:
                row[f"datalink_{c}"] = d.get(c)

            rows.append(row)

    except Exception as e:
        print(f"  ERROR: {e}")
        rows.append({
            "source_index": i,
            "target_name": target_name,
            "obs_id": obs_id,
            "obs_title": obs_title,
            "t_min": r.get("t_min"),
            "obscore_access_url": access_url,
            "error": repr(e),
        })

out = pd.DataFrame(rows)
out_csv = OUTDIR / "eso_espresso_hd128621_datalinks_expanded.csv"
out.to_csv(out_csv, index=False)

print("\nExpanded rows:", len(out))
print("Saved:", out_csv)

pd.set_option("display.max_colwidth", 180)
pd.set_option("display.max_rows", 200)

print("\nColumns:")
print(list(out.columns))

# Heurística para localizar URLs descargables
url_cols = [c for c in out.columns if "access" in c.lower() or "url" in c.lower()]
print("\nURL-like columns:", url_cols)

for c in url_cols:
    print(f"\nPreview URL column: {c}")
    print(out[c].dropna().astype(str).head(20).to_string(index=False))

# Filtrado probable de productos science FITS
mask = pd.Series(False, index=out.index)

for c in out.columns:
    lc = c.lower()
    if any(k in lc for k in ["semantics", "description", "content_type", "access_url", "url", "id"]):
        s = out[c].fillna("").astype(str)
        mask |= s.str.contains("fits|science|this|data|spectrum|SOBF|SFLX|ADP", case=False, regex=True)

cand = out[mask].copy()
cand_csv = OUTDIR / "eso_espresso_hd128621_download_candidates.csv"
cand.to_csv(cand_csv, index=False)

print("\nCandidate rows:", len(cand))
print("Saved:", cand_csv)

print("\nCandidate preview:")
show_cols = []
for c in [
    "target_name", "obs_title", "t_min", "t_exptime", "em_res_power",
    "datalink_semantics", "datalink_description", "datalink_content_type",
    "datalink_access_url", "datalink_access_format", "datalink_service_def"
]:
    if c in cand.columns:
        show_cols.append(c)

if len(cand) and show_cols:
    print(cand[show_cols].head(120).to_string(index=False))

print("\nDONE")
