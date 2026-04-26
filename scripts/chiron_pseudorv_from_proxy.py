from pathlib import Path
import pandas as pd
import numpy as np

C_MS = 299792458.0

inp = Path("data_interim/qc/chiron_dx_proxy_detrended.csv")
out = Path("data_interim/qc/chiron_pseudorv.csv")

df = pd.read_csv(inp).copy()

if "dx_proxy_detrended" not in df.columns:
    raise SystemExit("missing dx_proxy_detrended")
if "x_center_mass" not in df.columns:
    raise SystemExit("missing x_center_mass")

x_ref = df["x_center_mass"].astype(float).median()

df["pseudo_rv_ms"] = C_MS * (df["dx_proxy_detrended"].astype(float) / x_ref)

mad = np.median(np.abs(df["dx_proxy_detrended"] - np.median(df["dx_proxy_detrended"])))
sigma_dx = 1.4826 * mad
df["pseudo_rv_err_ms"] = C_MS * (sigma_dx / x_ref)

df["pseudo_rv_ms"] = df["pseudo_rv_ms"] - df["pseudo_rv_ms"].mean()

keep = [
    "stem","date_obs","mjd_obs","bjd_tdb",
    "airmass","flux_total","x_center_mass","y_center_mass",
    "dx_pix_proxy_demean","dx_proxy_detrended",
    "pseudo_rv_ms","pseudo_rv_err_ms",
    "quality_flag","notes","file_spectrum"
]
for c in keep:
    if c not in df.columns:
        df[c] = np.nan

df = df[keep].sort_values("date_obs")
df.to_csv(out, index=False)

print(df[["date_obs","dx_proxy_detrended","pseudo_rv_ms","pseudo_rv_err_ms"]].to_string(index=False))
print("\nRMS pseudo_rv_ms =", df["pseudo_rv_ms"].std())
print(f"\nSaved {out}")
