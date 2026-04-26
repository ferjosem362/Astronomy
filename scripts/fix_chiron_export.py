from pathlib import Path
import pandas as pd
import csv
from io import StringIO

inp = Path("data_raw/metadata/chiron_archive_export.csv")
out = Path("data_raw/metadata/chiron_archive_export_fixed.csv")

# Leemos con pandas normal
df = pd.read_csv(inp)

# Si ya estuviera bien, lo dejamos pasar
if df.iloc[:, 1:].notna().sum().sum() > 0:
    df.to_csv(out, index=False)
    print(f"already OK -> saved {out}")
    raise SystemExit

cols = list(df.columns)
fixed_rows = []

for _, row in df.iterrows():
    first = row.iloc[0]
    if pd.isna(first):
        continue
    parsed = next(csv.reader(StringIO(str(first))))
    # Si por alguna razón no cuadra, rellenamos/truncamos
    if len(parsed) < len(cols):
        parsed = parsed + [""] * (len(cols) - len(parsed))
    elif len(parsed) > len(cols):
        parsed = parsed[:len(cols)]
    fixed_rows.append(parsed)

fixed = pd.DataFrame(fixed_rows, columns=cols)
fixed.to_csv(out, index=False)
print(f"saved {out} rows={len(fixed)}")
print(fixed.head(10).to_string(index=False))
