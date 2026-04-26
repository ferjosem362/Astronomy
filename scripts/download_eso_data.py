from pathlib import Path
from astroquery.eso import Eso

OUTDIR = Path("data_raw/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)

eso = Eso()

targets = [
    "ALF CEN B",
    "ALPHA CEN B",
    "Alpha Cen B",
    "alpha cen b",
    "HD 128621",
    "HIP 71681",
]

for instrument in ["harps", "espresso"]:
    merged = None
    for target in targets:
        try:
            tbl = eso.query_instrument(
                instrument,
                column_filters={"target": f"= '{target}'"},
                authenticated=False,
            )
            if tbl is not None and len(tbl) > 0:
                if merged is None:
                    merged = tbl
                else:
                    try:
                        from astropy.table import vstack, unique
                        merged = unique(vstack([merged, tbl]), keys=list(tbl.colnames[:1]))
                    except Exception:
                        pass
                print(f"{instrument}: {target} -> {len(tbl)} rows")
        except Exception as e:
            print(f"{instrument}: {target} -> ERROR: {e}")

    if merged is None or len(merged) == 0:
        print(f"No rows found for {instrument}")
        continue

    out = OUTDIR / f"{instrument}_query_alpha_cen_b.ecsv"
    merged.write(out, format="ascii.ecsv", overwrite=True)
    print(f"Saved {len(merged)} rows to {out}")
