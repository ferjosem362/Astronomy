import csv
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
import socket
import time

socket.setdefaulttimeout(20)

inp = Path("data_raw/metadata/alpcenb_espresso_only.csv")
out_this = Path("data_raw/metadata/alpcenb_espresso_this_links.csv")
out_aux  = Path("data_raw/metadata/alpcenb_espresso_aux_links.csv")
out_err  = Path("data_raw/metadata/alpcenb_espresso_errors.csv")

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

with out_this.open("w", newline="") as f_this, \
     out_aux.open("w", newline="") as f_aux, \
     out_err.open("w", newline="") as f_err:

    w_this = csv.writer(f_this)
    w_aux  = csv.writer(f_aux)
    w_err  = csv.writer(f_err)

    header = [
        "target_name","instrument_name","dp_id","datalink_url",
        "file_url","content_type","content_length","eso_origfile","eso_category","semantics"
    ]
    w_this.writerow(header)
    w_aux.writerow(header)
    w_err.writerow(["dp_id","access_url","error"])

    for i, row in enumerate(base_rows, start=1):
        dp_id = row["dp_id"]
        print(f"[{i}/{len(base_rows)}] {dp_id}", flush=True)
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
            f_this.flush()
            f_aux.flush()
            print(f"  -> this={n_this}, aux={n_aux}", flush=True)
        except Exception as e:
            w_err.writerow([dp_id, row["access_url"], str(e)])
            f_err.flush()
            print(f"  ERROR: {e}", flush=True)
        time.sleep(0.2)

print("saved", out_this)
print("saved", out_aux)
print("saved", out_err)
