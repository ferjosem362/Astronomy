import csv
import subprocess
from pathlib import Path

downloads = [
    ("data_raw/metadata/alpcenb_harps_this_links.csv", "data_raw/harps"),
    ("data_raw/metadata/alpcenb_harps_aux_links.csv", "data_raw/harps_aux"),
]

for csv_path, outdir in downloads:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(csv_path) as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        url = row["file_url"].strip()
        if not url.startswith("http"):
            continue
        name = row["eso_origfile"].strip() or (row["dp_id"].strip() + ".dat")
        outfile = outdir / name
        print("Downloading", outfile.name)
        subprocess.run(["wget", "-nc", "-O", str(outfile), url], check=False)
