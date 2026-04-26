#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

INV = Path("data_work/serval/espresso_s1d/serval_inventory.csv")
OUT = Path("data_work/serval/espresso_s1d")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INV)

# Normalizar nombres por si el inventario cambió algo.
def has_col(c):
    return c in df.columns

if not has_col("resolved_path"):
    raise SystemExit("No encuentro columna resolved_path en serval_inventory.csv")

text = df.astype(str).apply(lambda r: " ".join(r.values), axis=1).str.upper()

is_acb = (
    text.str.contains("ALFCENB", na=False)
    | text.str.contains("ALF CEN B", na=False)
    | text.str.contains("HD128621", na=False)
)

is_aca = (
    text.str.contains("ALFCEN_ES_SFLX", na=False)
    | text.str.contains("ALF CEN ", na=False) & ~is_acb
)

ins = df["ins_mode"].fillna("").astype(str).str.upper() if "ins_mode" in df.columns else pd.Series("", index=df.index)
pro = df["pro_catg"].fillna("").astype(str).str.upper() if "pro_catg" in df.columns else pd.Series("", index=df.index)
hint = df["product_hint"].fillna("").astype(str).str.upper() if "product_hint" in df.columns else pd.Series("", index=df.index)

exists = df["exists"].astype(bool) if "exists" in df.columns else df["resolved_path"].map(lambda p: Path(str(p)).exists())

is_hr = ins.str.contains("SINGLEHR|HR", regex=True, na=False) & ~ins.str.contains("UHR", na=False)
is_uhr = ins.str.contains("UHR", na=False)
is_s1d = pro.str.contains("S1D", na=False) | hint.str.contains("ES_SOBF|ES_SFLX", regex=True, na=False) | text.str.contains("S1D", na=False)

mjd = pd.to_numeric(df["mjd_obs"], errors="coerce") if "mjd_obs" in df.columns else pd.Series(np.nan, index=df.index)

# Listas útiles.
lists = {
    "filelist_acb_all_s1d.txt": exists & is_acb & is_s1d,
    "filelist_acb_hr_all.txt": exists & is_acb & is_s1d & is_hr,
    "filelist_acb_hr_2023plus.txt": exists & is_acb & is_s1d & is_hr & (mjd >= 60042),
    "filelist_acb_hr_2023plus_no_early_offset_block.txt": exists & is_acb & is_s1d & is_hr & (mjd >= 60086),
    "filelist_acb_uhr.txt": exists & is_acb & is_s1d & is_uhr,
}

summary = []

for name, mask in lists.items():
    paths = df.loc[mask, "resolved_path"].dropna().astype(str).tolist()
    path = OUT / name
    path.write_text("\n".join(paths) + ("\n" if paths else ""))
    summary.append({
        "filelist": name,
        "n": len(paths),
        "mjd_min": float(mjd[mask].min()) if len(paths) and np.isfinite(mjd[mask]).any() else np.nan,
        "mjd_max": float(mjd[mask].max()) if len(paths) and np.isfinite(mjd[mask]).any() else np.nan,
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUT / "serval_clean_filelists_summary.csv", index=False)

print(summary_df.to_string(index=False))
print()
print("Recommended first SERVAL list:")
print(OUT / "filelist_acb_hr_2023plus.txt")
print()
print("Conservative post-offset-block list:")
print(OUT / "filelist_acb_hr_2023plus_no_early_offset_block.txt")
print()
print("DONE")
