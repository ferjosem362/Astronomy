from pathlib import Path
import pandas as pd
import numpy as np

inv = pd.read_csv("data_interim/qc/chiron_master_bjd.csv").copy()
feat = pd.read_csv("data_interim/qc/chiron_raw_features.csv").copy()
dx = pd.read_csv("data_interim/qc/chiron_order_pseudorv_v2.csv").copy()

# asegurar stem en todos
for df in [inv, feat, dx]:
    if "stem" not in df.columns:
        if "file_spectrum" in df.columns:
            df["stem"] = (
                df["file_spectrum"].astype(str)
                .str.split("/").str[-1]
                .str.replace(".fits.fz", "", regex=False)
            )
        elif "file" in df.columns:
            df["stem"] = (
                df["file"].astype(str)
                .str.split("/").str[-1]
                .str.replace(".fits.fz", "", regex=False)
            )

m = inv.merge(
    feat,
    on="stem",
    how="left",
    suffixes=("_inv", "_feat")
).merge(
    dx[["stem","dx_combined","dx_scatter","n_orders_used","file"]],
    on="stem",
    how="left",
    suffixes=("", "_dx")
)

print("COLUMNS AFTER MERGE:")
print(list(m.columns))

out = pd.DataFrame()
out["instrument"] = "CHIRON"
out["stem"] = m["stem"]
out["target_name"] = m["target_name"] if "target_name" in m.columns else np.nan
out["date_obs"] = m["date_obs"] if "date_obs" in m.columns else np.nan
out["bjd_tdb"] = m["bjd_tdb"] if "bjd_tdb" in m.columns else np.nan
out["mjd_obs"] = m["mjd_obs"] if "mjd_obs" in m.columns else np.nan
out["rv_rel"] = m["dx_combined"]
out["rv_rel_err"] = m["dx_scatter"]
out["airmass"] = m["airmass"] if "airmass" in m.columns else np.nan
out["flux_total"] = m["flux_total"] if "flux_total" in m.columns else np.nan
out["median_counts"] = m["median_counts"] if "median_counts" in m.columns else np.nan
out["p95_counts"] = m["p95_counts"] if "p95_counts" in m.columns else np.nan
out["x_center_mass"] = m["x_center_mass"] if "x_center_mass" in m.columns else np.nan
out["y_center_mass"] = m["y_center_mass"] if "y_center_mass" in m.columns else np.nan
out["n_orders_used"] = m["n_orders_used"]

# ruta del espectro: prueba varias columnas
if "file_spectrum" in m.columns:
    out["file_spectrum"] = m["file_spectrum"]
elif "file" in m.columns:
    out["file_spectrum"] = m["file"]
elif "file_inv" in m.columns:
    out["file_spectrum"] = m["file_inv"]
elif "file_feat" in m.columns:
    out["file_spectrum"] = m["file_feat"]
else:
    out["file_spectrum"] = np.nan

out["notes"] = "CHIRON v2 relative pseudoRV from order-derivative matching"

out.to_csv("data_interim/qc/chiron_v2_master.csv", index=False)
print("\nPREVIEW:")
print(out.head(20).to_string(index=False))
print(f"\nSaved data_interim/qc/chiron_v2_master.csv rows={len(out)}")
