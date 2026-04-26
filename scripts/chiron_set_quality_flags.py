from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/chiron_master_terms.csv")
out = Path("data_interim/qc/chiron_master_terms_qc.csv")

df = pd.read_csv(inp).copy()

flag = np.zeros(len(df), dtype=np.int32)
notes = df["notes"].fillna("").astype(str)

# bit 0 -> saturación
if "sat_frac" in df.columns:
    sat = (df["sat_frac"].fillna(0).to_numpy() > 0.01).astype(np.int32)
    flag |= sat

# bit 1 -> airmass alta
if "airmass" in df.columns:
    high_airmass = (df["airmass"].fillna(0).to_numpy() > 1.30).astype(np.int32)
    flag |= (high_airmass << 1)

# bit 2 -> flujo anómalo
if "flux_total" in df.columns:
    flux = pd.to_numeric(df["flux_total"], errors="coerce")
    med = flux.median()
    bad_flux = ((flux < 0.8 * med) | (flux > 1.2 * med)).fillna(False).to_numpy().astype(np.int32)
    flag |= (bad_flux << 2)

df["quality_flag"] = flag
df["notes"] = notes

df.to_csv(out, index=False)

print(df[["stem", "date_obs", "airmass", "quality_flag"]].head(20).to_string(index=False))
print("\nFlag counts:")
print(df["quality_flag"].value_counts(dropna=False).sort_index().to_string())
print(f"\nSaved {out}")
