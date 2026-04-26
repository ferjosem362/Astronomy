from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/chiron_v2_master.csv")
out = Path("data_interim/qc/chiron_v2_decorrelated.csv")

df = pd.read_csv(inp).copy()
df = df[df["rv_rel"].notna() & df["bjd_tdb"].notna()].copy()

features = ["airmass","flux_total","median_counts","p95_counts","x_center_mass","y_center_mass"]
X = df[features].copy()

for c in features:
    X[c] = pd.to_numeric(X[c], errors="coerce")
    med = X[c].median()
    if pd.isna(med):
        med = 0.0
    X[c] = X[c].fillna(med)

use = []
for c in features:
    s = float(np.nanstd(X[c].values))
    if np.isfinite(s) and s > 0:
        use.append(c)

X = X[use].copy()

for c in use:
    mu = X[c].mean()
    sd = X[c].std()
    if np.isfinite(sd) and sd > 0:
        X[c] = (X[c] - mu) / sd
    else:
        X[c] = 0.0

y = pd.to_numeric(df["rv_rel"], errors="coerce").values
e = pd.to_numeric(df["rv_rel_err"], errors="coerce").values

bad = ~np.isfinite(e) | (e <= 0)
if bad.any():
    e[bad] = np.nanmedian(e[~bad]) if (~bad).any() else 1.0

w = 1.0 / (e**2)
sw = np.sqrt(w)

A = np.column_stack([np.ones(len(X))] + [X[c].values for c in X.columns])
beta, *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)

model = A @ beta
resid = y - model

df["rv_model_decor"] = model
df["rv_resid_decor"] = resid

out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print("Used features:", list(X.columns))
print("Beta:", beta)
print(df[["date_obs","rv_rel","rv_model_decor","rv_resid_decor"]].to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
