from pathlib import Path
import pandas as pd

master_fp = Path("data_interim/qc/multinstrument_master.csv")
proxy_fp = Path("data_interim/qc/chiron_dx_proxy_detrended.csv")
out_fp = Path("data_interim/qc/multinstrument_master_with_chiron_proxy.csv")

master = pd.read_csv(master_fp)
proxy = pd.read_csv(proxy_fp)

keep = proxy[["stem","dx_pix_proxy","dx_pix_proxy_demean","dx_proxy_detrended","ccf_peak_proxy"]].copy()
out = master.merge(keep, on="stem", how="left")

out.to_csv(out_fp, index=False)

print(out[out["instrument"]=="CHIRON"][[
    "instrument","stem","date_obs","dx_pix_proxy","dx_pix_proxy_demean","dx_proxy_detrended"
]].head(20).to_string(index=False))
print(f"\nSaved {out_fp}")
