from pathlib import Path
import pandas as pd
import numpy as np
from astropy.io import fits

inp = Path("data_interim/qc/chiron_master_bjd.csv")
out = Path("data_interim/qc/chiron_raw_features.csv")

df = pd.read_csv(inp).copy()
rows = []

for _, row in df.iterrows():
    fp = Path(row["file_spectrum"])
    rec = {
        "stem": row["stem"],
        "file_spectrum": str(fp),
    }
    try:
        with fits.open(fp) as hdul:
            img = hdul[1].data.astype(float)

        finite = np.isfinite(img)
        vals = img[finite]

        rec["median_counts"] = float(np.nanmedian(vals))
        rec["p95_counts"] = float(np.nanpercentile(vals, 95))
        rec["max_counts"] = float(np.nanmax(vals))
        rec["sat_frac"] = float(np.mean(vals > 60000))
        rec["flux_total"] = float(np.nansum(vals))

        ny, nx = img.shape
        y0, y1 = int(0.25*ny), int(0.75*ny)
        x0, x1 = int(0.25*nx), int(0.75*nx)
        center = img[y0:y1, x0:x1]
        rec["flux_center"] = float(np.nansum(center))

        # centro de masa simple
        yy, xx = np.indices(img.shape)
        w = np.where(np.isfinite(img), np.clip(img, 0, None), 0.0)
        wsum = w.sum()
        if wsum > 0:
            rec["x_center_mass"] = float((xx * w).sum() / wsum)
            rec["y_center_mass"] = float((yy * w).sum() / wsum)
        else:
            rec["x_center_mass"] = np.nan
            rec["y_center_mass"] = np.nan

    except Exception as e:
        rec["error"] = str(e)

    rows.append(rec)

feat = pd.DataFrame(rows)
feat.to_csv(out, index=False)
print(feat.head(10).to_string(index=False))
print(f"\nSaved {out}")
