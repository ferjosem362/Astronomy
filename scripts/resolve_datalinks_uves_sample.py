import csv
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
import socket

socket.setdefaulttimeout(20)

inp = Path("data_raw/metadata/alpcenb_uves_sample.csv")
out_this = Path("data_raw/metadata/alpcenb_uves_sample_this_links.csv")
out_aux  = Path("data_raw/metadata/alpcenb_uves_sample_aux_links.csv")

ns = {"v": "http://www.ivoa.net/xml/VOTable/v1.3"}

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

with out_this.open("w", newline="") as f_this, out_aux.open("w", newline="") as f_aux:
    w_this = csv.writer(f_this)
    w_aux  = csv.writer(f_aux)

    header = [
        "target_name","instrument_name","dp_id","datalink_url",
        "file_url","content_type","content_length","eso_origfile","eso_category","semantics"
    ]
    w_this.writerow(header)
    w_aux.writerow(header)

    for i, row in enumerate(base_rows, start=1):
        print(f"[{i}/{len(base_rows)}] {row['dp_id']}", flush=True)
        try:
            dl_rows = parse_datalink(row["access_url"])
            n_this = 0
            n_aux = 0
            for d in dl_rows:
                rec = [
                    row["target_name"], row["instrument_name"], row["dp_id"], row["access_url"],
                    d.get("access_url",""), d.get("content_type",""), d.get("content_length",""),
                    d.get("eso_origfile",""), d.get("eso_category",""), d.get("semantics","")
                ]
                if d.get("semantics") == "#this":
                    w_this.writerow(rec)
                    n_this += 1
                elif d.get("semantics") == "#auxiliary":
                    w_aux.writerow(rec)
                    n_aux += 1
            print(f"  -> this={n_this}, aux={n_aux}", flush=True)
        except Exception as e:
            print(f"  ERROR: {e}", flush=True)

print("saved", out_this)
print("saved", out_aux)
