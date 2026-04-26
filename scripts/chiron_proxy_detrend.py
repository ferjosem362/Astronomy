from pathlib import Path
import pandas as pd
import statsmodels.api as sm

inp = Path("data_interim/qc/chiron_dx_proxy.csv")
out = Path("data_interim/qc/chiron_dx_proxy_detrended.csv")

df = pd.read_csv(inp).copy()

use = ["airmass", "flux_total", "x_center_mass", "y_center_mass"]
use = [c for c in use if c in df.columns]

X = df[use].copy()
X = sm.add_constant(X, has_constant="add")
y = df["dx_pix_proxy_demean"].astype(float)

model = sm.OLS(y, X, missing="drop").fit()
pred = model.predict(X)
df["dx_proxy_detrended"] = y - pred

df.to_csv(out, index=False)

print(model.summary())
print()
print(df[["stem","date_obs","dx_pix_proxy_demean","dx_proxy_detrended"]].to_string(index=False))
print(f"\nSaved {out}")
