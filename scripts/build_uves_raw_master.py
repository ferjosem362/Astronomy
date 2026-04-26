from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_interim/qc/uves_inventory.csv")
out = Path("data_interim/qc/uves_master_raw.csv")

df = pd.read_csv(inp).copy()

txt = (
    df.fillna("")
      .astype(str)
      .agg(" ".join, axis=1)
      .str.lower()
)

mask_b = (
    txt.str.contains("alpha-cen-b", na=False) |
    txt.str.contains("alpha cen b", na=False) |
    txt.str.contains("alph cen b", na=False) |
    txt.str.contains("hd128621", na=False) |
    txt.str.contains("hd 128621", na=False)
)

df = df.loc[mask_b].copy().reset_index(drop=True)

def first_existing(frame, cols):
    out = pd.Series(np.nan, index=frame.index, dtype=object)
    for c in cols:
        if c in frame.columns:
            out = out.where(out.notna(), frame[c])
    return out

def tag_type(fn):
    fn = str(fn)
    if "UV_S2DX_" in fn:
        return "S2DX"
    if "UV_SSXP_" in fn:
        return "SSXP"
    if fn.startswith("ADP.") or "ADP." in fn:
        return "ADP"
    if "UV_SRED_" in fn:
        return "SRED"
    if "UV_SOBR_" in fn:
        return "SOBR"
    return "OTHER"

prod_tag = df["filename"].astype(str).map(tag_type)

target_name = first_existing(df, ["object0", "object1", "object2"])
date_obs    = first_existing(df, ["date_obs0", "date_obs1", "date_obs2"])
mjd_obs     = first_existing(df, ["mjd_obs0", "mjd_obs1", "mjd_obs2"])
airmass     = first_existing(df, ["airmass0", "airmass1", "airmass2"])
snr         = first_existing(df, ["snr0", "snr1", "snr2"])

outdf = pd.DataFrame(index=df.index)
outdf["instrument"] = "UVES"
outdf["stem"] = (
    df["filename"].astype(str)
      .str.replace(".fits", "", regex=False)
      .str.replace(".fz", "", regex=False)
)
outdf["target_name"] = target_name
outdf["date_obs"] = date_obs
outdf["bjd_tdb"] = np.nan
outdf["mjd_obs"] = pd.to_numeric(mjd_obs, errors="coerce")
outdf["rv"] = np.nan
outdf["rv_err"] = np.nan
outdf["berv"] = np.nan
outdf["fwhm"] = np.nan
outdf["contrast"] = np.nan
outdf["bis_span"] = np.nan
outdf["snr"] = pd.to_numeric(snr, errors="coerce")
outdf["airmass"] = pd.to_numeric(airmass, errors="coerce")
outdf["drift"] = np.nan
outdf["aca_sep_arcsec"] = np.nan
outdf["aca_contam_proxy"] = np.nan
outdf["file_spectrum"] = df["file"]
outdf["file_ccf"] = np.nan
outdf["file_bis"] = np.nan
outdf["file_e2ds"] = np.where(prod_tag == "S2DX", df["file"], np.nan)
outdf["file_drift"] = np.nan
outdf["quality_flag"] = 0
outdf["notes"] = "UVES raw only; RV not yet extracted; prod_tag=" + prod_tag.astype(str)

out.parent.mkdir(parents=True, exist_ok=True)
outdf.to_csv(out, index=False)

print(outdf[["instrument","stem","target_name","date_obs","mjd_obs","notes"]].head(10).to_string(index=False))
print(f"\nSaved {out} rows={len(outdf)}")
