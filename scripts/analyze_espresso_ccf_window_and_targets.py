#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

INP = Path("data_interim/qc/periodogram_inputs/espresso_eso_ccf_curated_hr_detrended.csv")
PG = Path("results/periodograms_espresso_eso_ccf_curated_hr_detrended/gls_periodogram.csv")
OUTDIR = Path("results/periodograms_espresso_eso_ccf_curated_hr_detrended")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INP)
df["time"] = pd.to_numeric(df["time"], errors="coerce")
df["rv"] = pd.to_numeric(df["rv"], errors="coerce")
df["rv_err"] = pd.to_numeric(df["rv_err"], errors="coerce")
df = df[np.isfinite(df["time"]) & np.isfinite(df["rv"]) & np.isfinite(df["rv_err"])].copy()

t = df["time"].to_numpy()
y = df["rv"].to_numpy()
dy = df["rv_err"].to_numpy()

pg = pd.read_csv(PG)
pg["period_d"] = pd.to_numeric(pg["period_d"], errors="coerce")
pg["frequency_1_d"] = pd.to_numeric(pg["frequency_1_d"], errors="coerce")
pg["power"] = pd.to_numeric(pg["power"], errors="coerce")
pg = pg[np.isfinite(pg["period_d"]) & np.isfinite(pg["frequency_1_d"]) & np.isfinite(pg["power"])].copy()

freq = pg["frequency_1_d"].to_numpy()
period = pg["period_d"].to_numpy()

# Window function: constant signal sampled at t
w_y = np.ones_like(t)
w_ls = LombScargle(t, w_y, fit_mean=False, center_data=False)
freq = np.asarray(freq, dtype=float)
freq = freq[np.isfinite(freq)]
freq = np.unique(freq)
freq = np.sort(freq)

if len(freq) < 3:
    raise RuntimeError(f"Frequency grid demasiado corta: {len(freq)}")

if not np.all(np.diff(freq) > 0):
    raise RuntimeError("Frequency grid no es estrictamente creciente")

w_power = w_ls.power(freq)
win = pd.DataFrame({"period_d": period, "frequency_1_d": freq, "window_power": w_power})
win.to_csv(OUTDIR / "window_periodogram.csv", index=False)

# top window peaks
p = win["window_power"].to_numpy()
loc = np.where((p[1:-1] > p[:-2]) & (p[1:-1] > p[2:]))[0] + 1
wpeaks = win.iloc[loc].sort_values("window_power", ascending=False).head(50).copy()
wpeaks.to_csv(OUTDIR / "window_top_peaks.csv", index=False)

# Targeted powers at exact periods + local best around them
targets = [
    0.5, 0.970108077359177, 0.9941836673034685, 1.0, 1.008383808261102,
    1.028872682954664, 2.0, 2.024, 2.257, 3.236, 5.0, 7.0,
    16.3, 18.3, 18.98, 26.66, 32.33, 35.67, 46.86, 55.27,
    89.44, 123.39, 170.44
]

ls = LombScargle(t, y, dy, fit_mean=True, center_data=True)

rows = []
for P0 in targets:
    f0 = 1.0 / P0
    pow0 = float(ls.power(f0))
    win0 = float(w_ls.power(f0))

    lo = max(0.5, P0 * 0.95)
    hi = min(454.31, P0 * 1.05)
    fgrid = np.linspace(1 / hi, 1 / lo, 20000)
    pgrid = ls.power(fgrid)
    ib = int(np.nanargmax(pgrid))
    bestP = float(1 / fgrid[ib])
    bestPow = float(pgrid[ib])
    bestWin = float(w_ls.power(fgrid[ib]))

    rows.append({
        "target_period_d": P0,
        "exact_power": pow0,
        "exact_window_power": win0,
        "local_best_period_d": bestP,
        "local_best_power": bestPow,
        "local_best_window_power": bestWin,
        "local_window_ratio": bestPow / bestWin if bestWin > 0 else np.nan,
    })

targets_df = pd.DataFrame(rows)
targets_df.to_csv(OUTDIR / "target_period_diagnostics.csv", index=False)

print("Saved:", OUTDIR / "window_periodogram.csv")
print("Saved:", OUTDIR / "window_top_peaks.csv")
print("Saved:", OUTDIR / "target_period_diagnostics.csv")

print("\nTop window peaks:")
print(wpeaks.head(25).to_string(index=False))

print("\nTarget diagnostics:")
print(targets_df.to_string(index=False))

print("\nDONE")
