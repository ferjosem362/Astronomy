#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

INP = Path("data_interim/qc/periodogram_inputs/espresso_eso_ccf_curated_hr_detrended.csv")
OUTDIR = Path("results/periodograms_espresso_eso_ccf_curated_hr_detrended")
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INP)

for c in ["time", "rv", "rv_err"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df[np.isfinite(df["time"]) & np.isfinite(df["rv"]) & np.isfinite(df["rv_err"])].copy()
df = df.sort_values("time").reset_index(drop=True)

t = df["time"].to_numpy()
y = df["rv"].to_numpy()
dy = df["rv_err"].to_numpy()

# Centrado final global, por seguridad
y = y - np.nanmedian(y)

min_period = 0.5
max_period = min(500.0, 0.9 * (np.nanmax(t) - np.nanmin(t)))

freq = np.linspace(1 / max_period, 1 / min_period, 200000)
ls = LombScargle(t, y, dy, fit_mean=True, center_data=True)
power = ls.power(freq)
period = 1 / freq

faps = {
    "fap_10": ls.false_alarm_level(0.10, minimum_frequency=freq.min(), maximum_frequency=freq.max()),
    "fap_01": ls.false_alarm_level(0.01, minimum_frequency=freq.min(), maximum_frequency=freq.max()),
    "fap_001": ls.false_alarm_level(0.001, minimum_frequency=freq.min(), maximum_frequency=freq.max()),
}

pg = pd.DataFrame({"period_d": period, "frequency_1_d": freq, "power": power})
pg = pg.sort_values("period_d").reset_index(drop=True)
pg.to_csv(OUTDIR / "gls_periodogram.csv", index=False)

# Peak finding simple: local maxima separados en periodo
p = pg["power"].to_numpy()
per = pg["period_d"].to_numpy()
loc = np.where((p[1:-1] > p[:-2]) & (p[1:-1] > p[2:]))[0] + 1

peaks = pg.iloc[loc].sort_values("power", ascending=False).copy()

selected = []
for _, r in peaks.iterrows():
    P = float(r["period_d"])
    if all(abs(P - s["period_d"]) / P > 0.01 for s in selected):
        selected.append(r.to_dict())
    if len(selected) >= 50:
        break

top = pd.DataFrame(selected)
top["fap_10_level"] = faps["fap_10"]
top["fap_01_level"] = faps["fap_01"]
top["fap_001_level"] = faps["fap_001"]
top.to_csv(OUTDIR / "gls_top_peaks.csv", index=False)

summary = {
    "input": str(INP),
    "n": len(df),
    "time_min": float(np.nanmin(t)),
    "time_max": float(np.nanmax(t)),
    "baseline_d": float(np.nanmax(t) - np.nanmin(t)),
    "rv_std_mps": float(np.nanstd(y, ddof=1)),
    "rv_err_median_mps": float(np.nanmedian(dy)),
    "min_period_d": min_period,
    "max_period_d": max_period,
    **faps,
}

pd.Series(summary).to_csv(OUTDIR / "summary.csv")

print("Input:", INP)
print("Rows:", len(df))
print("Baseline d:", summary["baseline_d"])
print("RV std m/s:", summary["rv_std_mps"])
print("RV err median m/s:", summary["rv_err_median_mps"])
print("Saved:", OUTDIR / "gls_periodogram.csv")
print("Saved:", OUTDIR / "gls_top_peaks.csv")
print("Saved:", OUTDIR / "summary.csv")
print()
print("FAP levels:", faps)
print()
print("Top peaks:")
print(top.head(20).to_string(index=False))

# Plot periodogram
plt.figure(figsize=(12, 6))
plt.plot(pg["period_d"], pg["power"], lw=0.8)
for name, val in faps.items():
    plt.axhline(val, ls="--", lw=0.8, label=name)
plt.xscale("log")
plt.xlabel("Period [days]")
plt.ylabel("GLS power")
plt.title("ESPRESSO ESO CCF curated HR detrended GLS")
plt.legend()
plt.tight_layout()
plt.savefig(OUTDIR / "gls_periodogram.png", dpi=180)
plt.close()

# Plot time series
plt.figure(figsize=(12, 5))
for inst, g in df.groupby("model_inst"):
    plt.errorbar(g["time"], g["rv"], yerr=g["rv_err"], fmt=".", ms=3, alpha=0.75, label=inst)
plt.xlabel("BJD")
plt.ylabel("Detrended RV [m/s]")
plt.title("ESPRESSO ESO CCF curated HR detrended RV")
plt.legend()
plt.tight_layout()
plt.savefig(OUTDIR / "timeseries_detrended.png", dpi=180)
plt.close()

print("Saved:", OUTDIR / "gls_periodogram.png")
print("Saved:", OUTDIR / "timeseries_detrended.png")
print("\nDONE")
