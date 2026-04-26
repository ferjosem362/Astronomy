import csv
import subprocess
from pathlib import Path

inp = Path("data_raw/metadata/alpcenb_espresso_core_products_full.csv")
outdir = Path("data_raw/espresso")
outdir.mkdir(parents=True, exist_ok=True)

with inp.open() as f:
    rows = list(csv.DictReader(f))

for i, row in enumerate(rows, start=1):
    url = row["file_url"].strip()
    name = row["eso_origfile"].strip()
    if not url.startswith("http"):
        continue
    outfile = outdir / name
    print(f"[{i}/{len(rows)}] {name}", flush=True)
    subprocess.run(["wget", "-nc", "-O", str(outfile), url], check=False)
