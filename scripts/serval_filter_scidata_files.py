from pathlib import Path
import pandas as pd
from astropy.io import fits

IN = Path("data_work/serval/espresso_s1d/filelist_acb_hr_2023plus_no_early_offset_block.txt")
OUTDIR = Path("data_work/serval/espresso_s1d")
OUT = OUTDIR / "filelist_acb_hr_2023plus_no_early_SCIDATA.txt"
INV = OUTDIR / "serval_scidata_inventory.csv"

rows = []
good = []

for line in IN.read_text().splitlines():
    f = Path(line.strip())
    if not str(f):
        continue

    try:
        with fits.open(f, memmap=False) as hdul:
            extnames = []
            shapes = []
            for hdu in hdul:
                ext = hdu.header.get("EXTNAME", "")
                extnames.append(str(ext))
                try:
                    shapes.append(str(hdu.data.shape) if hdu.data is not None else "")
                except Exception:
                    shapes.append("")
            has_scidata = "SCIDATA" in extnames
            has_errdata = "ERRDATA" in extnames
            has_qualdata = "QUALDATA" in extnames
            object_name = hdul[0].header.get("OBJECT", "")
            date_obs = hdul[0].header.get("DATE-OBS", "")
            pro_catg = hdul[0].header.get("HIERARCH ESO PRO CATG", "")
            tpl_id = hdul[0].header.get("HIERARCH ESO TPL ID", "")
    except Exception as e:
        extnames = []
        shapes = []
        has_scidata = False
        has_errdata = False
        has_qualdata = False
        object_name = ""
        date_obs = ""
        pro_catg = ""
        tpl_id = ""
        rows.append({
            "path": str(f),
            "exists": f.exists(),
            "ok_open": False,
            "error": repr(e),
            "has_scidata": False,
            "has_errdata": False,
            "has_qualdata": False,
            "object": object_name,
            "date_obs": date_obs,
            "pro_catg": pro_catg,
            "tpl_id": tpl_id,
            "extnames": "",
            "shapes": "",
        })
        continue

    if has_scidata and has_errdata:
        good.append(str(f))

    rows.append({
        "path": str(f),
        "exists": f.exists(),
        "ok_open": True,
        "error": "",
        "has_scidata": has_scidata,
        "has_errdata": has_errdata,
        "has_qualdata": has_qualdata,
        "object": object_name,
        "date_obs": date_obs,
        "pro_catg": pro_catg,
        "tpl_id": tpl_id,
        "extnames": "|".join(extnames),
        "shapes": "|".join(shapes),
    })

df = pd.DataFrame(rows)
df.to_csv(INV, index=False)

OUT.write_text("\n".join(good) + ("\n" if good else ""))

print("Input:", IN)
print("Input n:", len(rows))
print("SCIDATA+ERRDATA n:", len(good))
print("Saved inventory:", INV)
print("Saved filelist:", OUT)
print()
print("By extnames:")
print(df.groupby(["has_scidata", "has_errdata", "has_qualdata"], dropna=False).size().reset_index(name="n").to_string(index=False))
print()
print("First non-SCIDATA:")
bad = df[~df["has_scidata"]].head(20)
print(bad[["path", "pro_catg", "extnames"]].to_string(index=False))
print()
print("DONE")
