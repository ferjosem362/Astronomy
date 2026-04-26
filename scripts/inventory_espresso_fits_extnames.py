from pathlib import Path
import pandas as pd
from astropy.io import fits

ROOTS = [
    Path("data_raw/espresso"),
    Path("data_work/serval/espresso_s1d/fits"),
]

rows = []

for root in ROOTS:
    if not root.exists():
        continue
    for f in sorted(root.rglob("*.fits")):
        try:
            with fits.open(f, memmap=False) as hdul:
                h0 = hdul[0].header
                extnames = []
                shapes = []
                for hdu in hdul:
                    extnames.append(str(hdu.header.get("EXTNAME", "")))
                    try:
                        shapes.append(str(hdu.data.shape) if hdu.data is not None else "")
                    except Exception:
                        shapes.append("")
                rows.append({
                    "path": str(f),
                    "ok": True,
                    "object": h0.get("OBJECT", ""),
                    "date_obs": h0.get("DATE-OBS", ""),
                    "pro_catg": h0.get("HIERARCH ESO PRO CATG", ""),
                    "tpl_id": h0.get("HIERARCH ESO TPL ID", ""),
                    "ins_mode": h0.get("HIERARCH ESO INS MODE", ""),
                    "has_scidata": "SCIDATA" in extnames,
                    "has_errdata": "ERRDATA" in extnames,
                    "has_qualdata": "QUALDATA" in extnames,
                    "has_spectrum": "SPECTRUM" in extnames,
                    "extnames": "|".join(extnames),
                    "shapes": "|".join(shapes),
                    "size_mb": f.stat().st_size / 1024 / 1024,
                })
        except Exception as e:
            rows.append({
                "path": str(f),
                "ok": False,
                "error": repr(e),
            })

df = pd.DataFrame(rows)
out = Path("data_raw/metadata/espresso_local_fits_ext_inventory.csv")
out.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out, index=False)

print("Total FITS:", len(df))
print("Saved:", out)
print()

if len(df):
    print("By pro_catg/ext type:")
    cols = ["pro_catg", "has_scidata", "has_spectrum"]
    print(df.groupby(cols, dropna=False).size().reset_index(name="n").sort_values("n", ascending=False).to_string(index=False))
    print()
    print("SCIDATA candidates:")
    cand = df[df.get("has_scidata", False).fillna(False)]
    print(cand[["path", "pro_catg", "date_obs", "extnames", "shapes"]].head(50).to_string(index=False))
    print()
    print("S2D-ish by name/pro_catg:")
    text = (
        df.get("path", "").astype(str).str.upper() + " " +
        df.get("pro_catg", "").astype(str).str.upper() + " " +
        df.get("extnames", "").astype(str).str.upper()
    )
    s2d = df[text.str.contains("S2D|SCIDATA|ERRDATA|QUALDATA", regex=True, na=False)]
    print(s2d[["path", "pro_catg", "date_obs", "extnames", "shapes"]].head(100).to_string(index=False))

print()
print("DONE")
