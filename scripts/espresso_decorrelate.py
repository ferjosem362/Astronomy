from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/multinstrument_master.csv")
out = Path("data_interim/qc/espresso_decorrelated.csv")

df = pd.read_csv(inp).copy()
df = df[df["instrument"] == "ESPRESSO"].copy()
df = df[df["rv"].notna() & df["bjd_tdb"].notna()].copy()

features = [
    "berv","fwhm","contrast","snr","drift",
    "year_sin_1","year_cos_1","year_sin_2","year_cos_2"
]

use0 = [c for c in features if c in df.columns]

X = df[use0].copy()
for c in use0:
    X[c] = pd.to_numeric(X[c], errors="coerce")
    med = X[c].median()
    if pd.isna(med):
        med = 0.0
    X[c] = X[c].fillna(med)

# quitar columnas sin varianza
use = []
dropped_constant = []
for c in use0:
    s = float(np.nanstd(X[c].values))
    if np.isfinite(s) and s > 0:
        use.append(c)
    else:
        dropped_constant.append(c)

X = X[use].copy()

# estandarizar
for c in use:
    mu = float(np.nanmean(X[c].values))
    s = float(np.nanstd(X[c].values))
    if np.isfinite(s) and s > 0:
        X[c] = (X[c] - mu) / s
    else:
        X[c] = 0.0

# quitar columnas casi duplicadas/casi colineales por correlación extrema
keep = []
dropped_collinear = []
corr_thresh = 0.9999

for c in X.columns:
    if not keep:
        keep.append(c)
        continue
    bad = False
    for k in keep:
        a = X[c].values
        b = X[k].values
        if np.nanstd(a) == 0 or np.nanstd(b) == 0:
            continue
        r = np.corrcoef(a, b)[0,1]
        if np.isfinite(r) and abs(r) >= corr_thresh:
            bad = True
            dropped_collinear.append(c)
            break
    if not bad:
        keep.append(c)

X = X[keep].copy()

y = df["rv"].astype(float).values
e = pd.to_numeric(df["rv_err"], errors="coerce").values
bad_e = ~np.isfinite(e) | (e <= 0)
if bad_e.any():
    med_e = np.nanmedian(e[~bad_e]) if (~bad_e).any() else 1.0
    e[bad_e] = med_e

w = 1.0 / np.maximum(e, 1e-6)**2
sw = np.sqrt(w)

A = np.column_stack([np.ones(len(X))] + [X[c].values for c in X.columns])
Aw = A * sw[:, None]
yw = y * sw

# resolución robusta por least squares
beta, resid_sum, rank, svals = np.linalg.lstsq(Aw, yw, rcond=None)

model = A @ beta
resid = y - model

df["rv_model_decor"] = model
df["rv_resid_decor"] = resid

out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print("Initial features:", use0)
print("Dropped constant:", dropped_constant)
print("Dropped collinear:", dropped_collinear)
print("Used features:", list(X.columns))
print("Rank:", rank)
print("Singular values:", svals)
print("Beta:", beta)

print("\nPreview:")
print(df[["date_obs","rv","rv_model_decor","rv_resid_decor"]].head(15).to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
