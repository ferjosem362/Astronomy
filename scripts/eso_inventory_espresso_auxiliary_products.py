#!/usr/bin/env python3
from pathlib import Path
import re
import pandas as pd

INCSV = Path("data_raw/metadata/eso_espresso_hd128621_datalinks_expanded.csv")
OUTDIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INCSV)

aux = df[df["datalink_semantics"].astype(str).eq("#auxiliary")].copy()
aux["size_mb"] = aux["datalink_content_length"].fillna(0) / 1024**2

def extract_id(url):
    m = re.search(r"((?:ADP|ESPRE)\.[A-Za-z0-9T:\.\-]+)", str(url))
    return m.group(1) if m else ""

aux["file_id"] = aux["datalink_access_url"].map(extract_id)

def token_flags(row):
    txt = " ".join(str(row.get(c, "")) for c in [
        "obs_title", "datalink_access_url", "datalink_eso_origfile",
        "datalink_eso_category", "datalink_description",
        "datalink_content_type"
    ]).upper()
    keys = ["CCF", "RV", "RVC", "BIS", "FWHM", "S2D", "S1D", "ORDER", "TELL", "BLAZE", "WAVE", "CALIB", "DRS", "QC"]
    return {f"has_{k.lower()}": (k in txt) for k in keys}

flags = pd.DataFrame([token_flags(r) for _, r in aux.iterrows()], index=aux.index)
aux = pd.concat([aux, flags], axis=1)

aux_out = OUTDIR / "eso_espresso_auxiliary_inventory.csv"
aux.to_csv(aux_out, index=False)

summary = (
    aux.groupby([
        "datalink_content_type",
        "datalink_eso_category",
        "datalink_content_length",
    ], dropna=False)
    .agg(
        n=("datalink_access_url", "count"),
        n_obs=("obs_id", "nunique"),
        first_obs_title=("obs_title", "first"),
        first_origfile=("datalink_eso_origfile", "first"),
        first_url=("datalink_access_url", "first"),
    )
    .reset_index()
    .sort_values(["n", "datalink_content_length"], ascending=[False, False])
)

summary_out = OUTDIR / "eso_espresso_auxiliary_inventory_summary.csv"
summary.to_csv(summary_out, index=False)

interesting_mask = (
    aux["has_ccf"] | aux["has_rv"] | aux["has_rvc"] | aux["has_bis"] |
    aux["has_fwhm"] | aux["has_s2d"] | aux["has_order"] | aux["has_qc"]
)
interesting = aux[interesting_mask].copy()
interesting_out = OUTDIR / "eso_espresso_auxiliary_interesting_candidates.csv"
interesting.to_csv(interesting_out, index=False)

print("Input rows:", len(df))
print("Auxiliary rows:", len(aux))
print("Saved:", aux_out)
print("Saved:", summary_out)
print("Saved:", interesting_out)
print()

print("Auxiliary content/category/size summary:")
pd.set_option("display.max_colwidth", 160)
print(summary.head(80).to_string(index=False))

print("\nInteresting token counts:")
for c in [c for c in aux.columns if c.startswith("has_")]:
    n = int(aux[c].sum())
    if n:
        print(f"{c}: {n}")

print("\nInteresting preview:")
cols = [
    "target_name", "obs_id", "obs_title", "t_min",
    "datalink_content_type", "datalink_content_length", "size_mb",
    "datalink_eso_category", "datalink_eso_origfile",
    "datalink_access_url",
]
cols = [c for c in cols if c in interesting.columns]
if len(interesting):
    print(interesting[cols].head(120).to_string(index=False))
else:
    print("No interesting candidates by token search.")

print("\nDONE")
