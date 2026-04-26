from pathlib import Path
import pandas as pd
import numpy as np

inp = Path("data_raw/metadata/chiron_archive_export_fixed.csv")
out = Path("data_interim/qc/chiron_manifest.csv")

df = pd.read_csv(inp)

outdf = pd.DataFrame()
outdf["instrument"] = "CHIRON"
outdf["target_name"] = df["object"] if "object" in df.columns else np.nan
outdf["date_obs"] = df["date_obs"] if "date_obs" in df.columns else np.nan
outdf["mjd_obs"] = np.nan
outdf["snr"] = np.nan
outdf["download_url"] = np.nan
outdf["filename"] = df["archive_filename"] if "archive_filename" in df.columns else np.nan
outdf["original_filename"] = df["original_filename"] if "original_filename" in df.columns else np.nan
outdf["filesize"] = df["filesize"] if "filesize" in df.columns else np.nan
outdf["ra"] = df["ra"] if "ra" in df.columns else np.nan
outdf["dec"] = df["dec"] if "dec" in df.columns else np.nan
outdf["exposure"] = df["exposure"] if "exposure" in df.columns else np.nan
outdf["obs_type"] = df["obs_type"] if "obs_type" in df.columns else np.nan
outdf["proc_type"] = df["proc_type"] if "proc_type" in df.columns else np.nan
outdf["prop_id"] = df["prop_id"] if "prop_id" in df.columns else np.nan
outdf["telescope"] = df["telescope"] if "telescope" in df.columns else np.nan

out.parent.mkdir(parents=True, exist_ok=True)
outdf.to_csv(out, index=False)

print(f"saved {out} rows={len(outdf)}")
print(outdf.head(20).to_string(index=False))
