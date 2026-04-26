#!/usr/bin/env python3
from pathlib import Path
import shutil
import pandas as pd
import numpy as np

INDEX = Path("data_raw/metadata/espresso_s1d_index.csv")
OUT = Path("data_work/serval/espresso_s1d")
FITS_OUT = OUT / "fits"
LOGDIR = Path("logs/serval")

OUT.mkdir(parents=True, exist_ok=True)
FITS_OUT.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

print("Using index:", INDEX)
df = pd.read_csv(INDEX)
print("Index shape:", df.shape)
print("Columns:", list(df.columns))

# Detect file/path column
candidate_path_cols = [
    "file", "local_path", "path", "filename", "fits_path"
]
path_col = None
for c in candidate_path_cols:
    if c in df.columns:
        vals = df[c].dropna().astype(str)
        if len(vals) and vals.str.contains(".fits", case=False, regex=False).any():
            path_col = c
            break

if path_col is None:
    raise SystemExit("No encuentro columna de path FITS en espresso_s1d_index.csv")

print("Detected path column:", path_col)

d = df.copy()

# Texto robusto para filtrar producto/instrumento
text = d.apply(
    lambda r: " ".join("" if pd.isna(x) else str(x) for x in r.values),
    axis=1
).str.upper()

# Preferimos S1D de Alpha Cen B / HD128621, no CCF, no PNG, no tar.
mask = text.str.contains("S1D", na=False) | text.str.contains("SPECTRUM", na=False)
mask &= ~text.str.contains("CCF", na=False)
mask &= ~text.str.contains(".TAR", na=False)
mask &= ~text.str.contains(".PNG", na=False)

# Si el filtro fuese demasiado agresivo, usar todo lo que parezca FITS.
sel = d[mask].copy()
if len(sel) == 0:
    print("WARNING: filtro S1D vacío; uso todos los FITS del índice.")
    sel = d.copy()

# Normalizar paths
def resolve_path(x):
    p = Path(str(x))
    if p.exists():
        return p
    # Si la columna solo trae filename, buscar en data_raw/espresso/this
    p2 = Path("data_raw/espresso/this") / p.name
    if p2.exists():
        return p2
    return p

sel["_src_path"] = sel[path_col].map(resolve_path)
sel["_exists"] = sel["_src_path"].map(lambda p: p.exists())

print("Selected rows before exists:", len(sel))
print("Existing FITS:", int(sel["_exists"].sum()))

missing = sel[~sel["_exists"]]
if len(missing):
    print("\nMissing first 20:")
    print(missing[[path_col, "_src_path"]].head(20).to_string(index=False))

sel = sel[sel["_exists"]].copy()

if len(sel) == 0:
    raise SystemExit("No hay FITS existentes para preparar SERVAL.")

# Ordenar por fecha si existe
for c in ["mjd_obs", "MJD-OBS", "date_obs", "DATE-OBS", "filename"]:
    if c in sel.columns:
        sel = sel.sort_values(c)
        break

# Copiar/symlink. Symlink ahorra espacio.
filelist = []
manifest_rows = []

for i, row in sel.iterrows():
    src = Path(row["_src_path"])
    dst = FITS_OUT / src.name

    if not dst.exists():
        try:
            dst.symlink_to(src.resolve())
        except Exception:
            shutil.copy2(src, dst)

    filelist.append(str(dst.resolve()))
    manifest_rows.append({
        "src": str(src),
        "dst": str(dst),
        "filename": src.name,
    })

filelist_path = OUT / "filelist.txt"
filelist_path.write_text("\n".join(filelist) + "\n")

manifest = pd.DataFrame(manifest_rows)
manifest_path = OUT / "manifest.csv"
manifest.to_csv(manifest_path, index=False)

print()
print("Prepared SERVAL workspace:")
print("  OUT:", OUT)
print("  FITS:", FITS_OUT)
print("  filelist:", filelist_path)
print("  manifest:", manifest_path)
print("  n_files:", len(filelist))

print("\nFirst files:")
for p in filelist[:10]:
    print(" ", p)

print("\nDONE")
