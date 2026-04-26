from pathlib import Path
import pandas as pd
import subprocess

inp = Path("data_interim/qc/chiron_manifest.csv")
outdir = Path("data_raw/chiron")
outdir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(inp)

for i, row in df.iterrows():
    url = str(row.get("download_url", "")).strip()
    name = str(row.get("filename", "")).strip()
    if not url.startswith("http") or not name:
        continue
    outfile = outdir / name
    print(f"[{i+1}/{len(df)}] {name}")
    subprocess.run(["wget", "-nc", "-O", str(outfile), url], check=False)
