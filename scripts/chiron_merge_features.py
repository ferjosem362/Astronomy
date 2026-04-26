from pathlib import Path
import pandas as pd

inp_master = Path("data_interim/qc/chiron_master_bjd.csv")
inp_feat = Path("data_interim/qc/chiron_raw_features.csv")
out = Path("data_interim/qc/chiron_master_enriched.csv")

m = pd.read_csv(inp_master)
f = pd.read_csv(inp_feat)

df = m.merge(f, on=["stem","file_spectrum"], how="left")

# proxy simple de calidad / contaminación instrumental
if "flux_total" in df.columns and "airmass" in df.columns:
    df["aca_contam_proxy"] = df["flux_total"] / df["flux_total"].median()

# SNR proxy provisional
if "p95_counts" in df.columns and "median_counts" in df.columns:
    df["snr"] = (df["p95_counts"] - df["median_counts"]).clip(lower=0)

df.to_csv(out, index=False)
print(df.head(10).to_string(index=False))
print(f"\nSaved {out}")
