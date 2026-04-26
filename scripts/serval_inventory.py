#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(".").resolve()
WORK = Path("data_work/serval/espresso_s1d")
WORK.mkdir(parents=True, exist_ok=True)

INDEX_CANDIDATES = [
    Path("data_raw/metadata/espresso_s1d_index.csv"),
    Path("data_raw/metadata/eso_espresso_hd128621_this_manifest.csv"),
]

idx_path = next((p for p in INDEX_CANDIDATES if p.exists()), None)
if idx_path is None:
    raise SystemExit("No encuentro espresso_s1d_index.csv ni manifest #this.")

df = pd.read_csv(idx_path)
print("Using index:", idx_path)
print("Index shape:", df.shape)
print("Columns:", list(df.columns))

def pick_col(candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

path_col = pick_col(["file", "local_path", "resolved_path"])
filename_col = pick_col(["filename", "local_filename"])
mjd_col = pick_col(["mjd_obs", "MJD-OBS", "MJD_OBS", "t_min"])
date_col = pick_col(["date_obs", "DATE-OBS", "DATE_OBS"])
ins_col = pick_col(["ins_mode", "HIERARCH ESO INS MODE", "hierarch_eso_ins_mode"])
pro_col = pick_col(["pro_catg", "HIERARCH ESO PRO CATG", "hierarch_eso_pro_catg", "PRODCATG"])
target_col = pick_col(["target", "target_name", "OBJECT", "object", "HIERARCH ESO OBS TARG NAME"])
hint_col = pick_col(["product_hint", "espresso_product_hint"])

if path_col is None:
    raise SystemExit("No encuentro columna de path local en el índice.")

out = pd.DataFrame()
out["source_index_row"] = np.arange(len(df))
out["source_index_file"] = str(idx_path)
out["raw_path"] = df[path_col].astype(str)

def resolve_path(p):
    pp = Path(str(p))
    if pp.is_absolute():
        return str(pp)
    return str((ROOT / pp).resolve())

out["resolved_path"] = out["raw_path"].map(resolve_path)
out["filename"] = (
    df[filename_col].astype(str)
    if filename_col else out["resolved_path"].map(lambda p: Path(p).name)
)

out["exists"] = out["resolved_path"].map(lambda p: Path(p).exists())
out["size_bytes"] = out["resolved_path"].map(lambda p: Path(p).stat().st_size if Path(p).exists() else np.nan)

out["mjd_obs"] = pd.to_numeric(df[mjd_col], errors="coerce") if mjd_col else np.nan
out["date_obs"] = df[date_col].astype(str) if date_col else ""
out["ins_mode"] = df[ins_col].astype(str) if ins_col else ""
out["pro_catg"] = df[pro_col].astype(str) if pro_col else ""
out["target"] = df[target_col].astype(str) if target_col else ""
out["product_hint"] = df[hint_col].astype(str) if hint_col else ""

# Inferencias robustas desde texto/fichero
text = out.astype(str).apply(lambda r: " ".join(r.values), axis=1).str.upper()

out["is_acb"] = (
    text.str.contains("ALFCENB", na=False)
    | text.str.contains("ALF CEN B", na=False)
    | text.str.contains("HD128621", na=False)
)
out["is_aca"] = (
    text.str.contains("ALFCEN_ES_", na=False)
    | (
        text.str.contains("ALF CEN", na=False)
        & ~out["is_acb"]
    )
)

out["is_uhr"] = text.str.contains("UHR", na=False)
out["is_hr"] = text.str.contains("HR", na=False) & ~out["is_uhr"]
out["is_s1d"] = (
    text.str.contains("S1D", na=False)
    | text.str.contains("ES_SOBF", na=False)
    | text.str.contains("ES_SFLX", na=False)
)

# Symlink workspace
fits_dir = WORK / "fits"
fits_dir.mkdir(parents=True, exist_ok=True)

linked = 0
for p in out.loc[out["exists"], "resolved_path"]:
    src = Path(p)
    dst = fits_dir / src.name
    if not dst.exists():
        try:
            dst.symlink_to(src)
            linked += 1
        except FileExistsError:
            pass

# Filelist bruto de todo lo existente
all_paths = out.loc[out["exists"], "resolved_path"].astype(str).tolist()
(WORK / "filelist.txt").write_text("\n".join(all_paths) + ("\n" if all_paths else ""))

out_path = WORK / "serval_inventory.csv"
out.to_csv(out_path, index=False)

print()
print("Saved:", out_path)
print("Saved:", WORK / "filelist.txt")
print("Symlinks created:", linked)
print()
print("Exists counts:")
print(out["exists"].value_counts(dropna=False).to_string())
print()
print("Target flags:")
print(out[["is_acb", "is_aca", "is_hr", "is_uhr", "is_s1d"]].sum().to_string())
print()
print("MJD range:", out["mjd_obs"].min(), out["mjd_obs"].max())
print()
print("Preview:")
print(out[["filename", "target", "product_hint", "ins_mode", "pro_catg", "mjd_obs", "exists", "is_acb", "is_aca", "is_hr", "is_uhr"]].head(20).to_string(index=False))
print()
print("DONE")
