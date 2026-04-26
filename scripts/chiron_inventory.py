from pathlib import Path
from astropy.io import fits
import pandas as pd

root = Path("data_raw/chiron")
rows = []

def get_any(hdr, keys):
    for k in keys:
        if k in hdr:
            return hdr[k]
    return None

for fp in sorted(root.glob("*.fits.fz")):
    try:
        with fits.open(fp) as hdul:
            hdr0 = hdul[0].header
            hdr1 = hdul[1].header if len(hdul) > 1 else {}

            rows.append({
                "file": str(fp),
                "filename": fp.name,
                "n_hdu": len(hdul),

                "object0": get_any(hdr0, ["OBJECT","OBJNAME","TARGNAME","TARGET"]),
                "object1": get_any(hdr1, ["OBJECT","OBJNAME","TARGNAME","TARGET"]),

                "date_obs0": get_any(hdr0, ["DATE-OBS","DATEOBS"]),
                "date_obs1": get_any(hdr1, ["DATE-OBS","DATEOBS"]),

                "mjd_obs0": get_any(hdr0, ["MJD-OBS","MJDOBS","JD","MJD"]),
                "mjd_obs1": get_any(hdr1, ["MJD-OBS","MJDOBS","JD","MJD"]),

                "exptime0": get_any(hdr0, ["EXPTIME","EXPOSURE"]),
                "exptime1": get_any(hdr1, ["EXPTIME","EXPOSURE"]),

                "ra0": get_any(hdr0, ["RA"]),
                "ra1": get_any(hdr1, ["RA"]),
                "dec0": get_any(hdr0, ["DEC"]),
                "dec1": get_any(hdr1, ["DEC"]),

                "airmass0": get_any(hdr0, ["AIRMASS"]),
                "airmass1": get_any(hdr1, ["AIRMASS"]),

                "instrume0": get_any(hdr0, ["INSTRUME"]),
                "instrume1": get_any(hdr1, ["INSTRUME"]),
                "telescop0": get_any(hdr0, ["TELESCOP"]),
                "telescop1": get_any(hdr1, ["TELESCOP"]),

                "shape1": None if hdul[1].data is None else str(hdul[1].data.shape),
            })
    except Exception as e:
        rows.append({
            "file": str(fp),
            "filename": fp.name,
            "error": str(e),
        })

df = pd.DataFrame(rows)
out = Path("data_interim/qc/chiron_inventory.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print(df.to_string(index=False))
print(f"\nSaved {out} rows={len(df)}")
