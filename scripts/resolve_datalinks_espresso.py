import csv
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

inp = Path("data_raw/metadata/alpcenb_espresso_only.csv")
out_this = Path("data_raw/metadata/alpcenb_espresso_this_links.csv")
out_aux  = Path("data_raw/metadata/alpcenb_espresso_aux_links.csv")

ns = {"v": "http://www.ivoa.net/xml/VOTable/v1.3"}

def parse_datalink(url):
    with urllib.request.urlopen(url, timeout=60) as r:
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

    for row in base_rows:
        try:
            dl_rows = parse_datalink(row["access_url"])
            for d in dl_rows:
                rec = [
                    row["target_name"], row["instrument_name"], row["dp_id"], row["access_url"],
                    d.get("access_url",""), d.get("content_type",""), d.get("content_length",""),
                    d.get("eso_origfile",""), d.get("eso_category",""), d.get("semantics","")
                ]
                if d.get("semantics") == "#this":
                    w_this.writerow(rec)
                elif d.get("semantics") == "#auxiliary":
                    w_aux.writerow(rec)
        except Exception as e:
            print("ERROR", row.get("dp_id"), e)

print("saved", out_this)
print("saved", out_aux)
