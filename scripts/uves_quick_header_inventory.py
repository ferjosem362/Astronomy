from pathlib import Path
from astropy.io import fits
import pandas as pd

inp = Path("data_interim/qc/uves_master_raw.csv")
out = Path("data_interim/qc/uves_quick_header_inventory.csv")

df = pd.read_csv(inp)
m = df["notes"].str.contains("prod_tag=SOBR|prod_tag=SRED", na=False)
df = df.loc[m].copy()

rows = []
for _, r in df.iterrows():
    fp = Path(r["file_spectrum"])
    rec = {
        "stem": r["stem"],
        "target_name": r["target_name"],
        "date_obs_csv": r["date_obs"],
        "mjd_obs_csv": r["mjd_obs"],
        "file_spectrum": str(fp),
    }
    try:
        with fits.open(fp) as hdul:
            rec["n_hdu"] = len(hdul)
            for i in range(min(len(hdul), 3)):
                hdr = hdul[i].header
                for k in [
                    "OBJECT","DATE-OBS","MJD-OBS","EXPTIME","RA","DEC","AIRMASS",
                    "SNR","SPEC_RES","BERV","RV","VCORR","HELCOR","BJD","JD"
                ]:
                    rec[f"hdu{i}_{k.replace('-','_')}"] = hdr.get(k)
    except Exception as e:
        rec["error"] = str(e)
    rows.append(rec)

outdf = pd.DataFrame(rows)
out.parent.mkdir(parents=True, exist_ok=True)
outdf.to_csv(out, index=False)

print(outdf.head(20).to_string(index=False))
print(f"\nSaved {out} rows={len(outdf)}")
