from pathlib import Path
from astropy.io import fits
import pandas as pd

rows = []

files = sorted(Path("data_raw/uves").rglob("*.fits*"))
for fp in files:
    rec = {
        "file": str(fp),
        "filename": fp.name,
    }
    try:
        with fits.open(fp) as hdul:
            rec["n_hdu"] = len(hdul)
            for i in range(min(len(hdul), 3)):
                hdr = hdul[i].header
                rec[f"object{i}"]   = hdr.get("OBJECT")
                rec[f"date_obs{i}"] = hdr.get("DATE-OBS")
                rec[f"mjd_obs{i}"]  = hdr.get("MJD-OBS")
                rec[f"exptime{i}"]  = hdr.get("EXPTIME")
                rec[f"ra{i}"]       = hdr.get("RA")
                rec[f"dec{i}"]      = hdr.get("DEC")
                rec[f"airmass{i}"]  = hdr.get("AIRMASS")
                rec[f"instrume{i}"] = hdr.get("INSTRUME")
                rec[f"telescop{i}"] = hdr.get("TELESCOP")
                rec[f"snr{i}"]      = hdr.get("SNR")
            if len(hdul) > 1 and hdul[1].data is not None:
                try:
                    rec["shape1"] = str(hdul[1].data.shape)
                except Exception:
                    pass
    except Exception as e:
        rec["error"] = str(e)

    rows.append(rec)

df = pd.DataFrame(rows)
out = Path("data_interim/qc/uves_inventory.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print(df.head(20).to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
