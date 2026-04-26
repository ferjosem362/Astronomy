import csv
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
import socket

socket.setdefaulttimeout(20)

inp = Path("data_raw/metadata/alpcenb_espresso_sample5.csv")
out = Path("data_raw/metadata/alpcenb_espresso_sample5_products.csv")

ns = {"v": "http://www.ivoa.net/xml/VOTable/v1.3"}

wanted_categories = {
    "SCIENCE.SPECTRUM",
    "ANCILLARY.CCF",
    "ANCILLARY.DRIFT",
}

def parse_datalink(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=20) as r:
        xml_data = r.read()
    root = ET.fromstring(xml_data)

    fields = [f.attrib.get("name", "") for f in root.findall(".//v:FIELD", ns)]
    rows = []
    for tr in root.findall(".//v:TR", ns):
        vals = [td.text if td.text is not None else "" for td in tr.findall("v:TD", ns)]
        rows.append(dict(zip(fields, vals)))
    return rows

with inp.open() as f:
    base_rows = list(csv.DictReader(f))

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "target_name","dp_id","t_min","snr",
        "file_url","content_type","content_length","eso_origfile","eso_category","semantics"
    ])

    for i, row in enumerate(base_rows, start=1):
        print(f"[{i}/{len(base_rows)}] {row['dp_id']}", flush=True)
        try:
            dl_rows = parse_datalink(row["access_url"])
            kept = 0
            for d in dl_rows:
                cat = d.get("eso_category","")
                if cat in wanted_categories:
                    w.writerow([
                        row["target_name"], row["dp_id"], row["t_min"], row["snr"],
                        d.get("access_url",""), d.get("content_type",""), d.get("content_length",""),
                        d.get("eso_origfile",""), cat, d.get("semantics","")
                    ])
                    kept += 1
            print(f"  -> kept {kept}", flush=True)
        except Exception as e:
            print(f"  ERROR: {e}", flush=True)

print(f"saved {out}")
