from pathlib import Path
from astropy.io import fits
import pandas as pd

ROOTS = [Path("data_raw/harps"), Path("data_raw/espresso")]
rows = []

for root in ROOTS:
    if not root.exists():
        continue
    for fp in root.rglob("*.fits"):
        try:
            hdr = fits.getheader(fp, 0)
            rows.append({
                "file": str(fp),
                "instrument": hdr.get("INSTRUME"),
                "object": hdr.get("OBJECT"),
                "date_obs": hdr.get("DATE-OBS"),
                "mjd_obs": hdr.get("MJD-OBS"),
                "exptime": hdr.get("EXPTIME"),
                "snr": hdr.get("SNR") or hdr.get("HIERARCH ESO DRS SPE EXT SN50"),
                "prog_id": hdr.get("HIERARCH ESO OBS PROG ID"),
                "ra": hdr.get("RA"),
                "dec": hdr.get("DEC"),
                "arcfile": hdr.get("ARCFILE"),
            })
        except Exception as e:
            rows.append({"file": str(fp), "error": str(e)})

df = pd.DataFrame(rows)
out = Path("docs/fits_inventory.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)
print(f"Saved {out} with {len(df)} rows")
