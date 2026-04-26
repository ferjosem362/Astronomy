#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

INPUT = "data_interim/qc/multinstrument_rvonly.csv"
OUTDIR = Path("results/periodograms_dashboard")
OUTDIR.mkdir(parents=True, exist_ok=True)

P_MIN = 0.5
P_MAX = 120.0
SAMPLES_PER_PEAK = 25
VETO_FRAC = 0.05

M10 = 11.867564152695252

MARKS = {
    "0.5 d": 0.5,
    "1 d": 1.0,
    "1.62 d": 1.62,
    "3.24 d": 3.24,
    "M10": M10,
    "P53": 53.4,
    "C70": 70.5,
    "P_rot 36.2": 36.2,
    "P_rot 38.7": 38.7,
    "P_rot/2": 36.2 / 2,
    "P_rot/2b": 38.7 / 2,
    "P_rot/3": 36.2 / 3,
    "P_rot/3b": 38.7 / 3,
    "P_rot/4": 36.2 / 4,
    "P_rot/4b": 38.7 / 4,
}

def classical_window(t, freq):
    tt = np.asarray(t, dtype=float) - np.nanmin(t)
    n = len(tt)
    out = np.empty_like(freq)
    chunk = 2000
    for i in range(0, len(freq), chunk):
        ff = freq[i:i + chunk]
        phase = np.exp(2j * np.pi * ff[:, None] * tt[None, :])
        out[i:i + chunk] = np.abs(np.sum(phase, axis=1)) ** 2 / n ** 2
    return out

def prepare(df):
    df = df.copy()
    if "quality_flag" in df.columns:
        df = df[df["quality_flag"].fillna(0).astype(int) == 0].copy()
    df = df.dropna(subset=["instrument", "bjd_tdb", "rv", "rv_err"]).copy()
    df["rv_centered"] = df["rv"] - df.groupby("instrument")["rv"].transform("median")
    return df

def calc(sub):
    t = sub["bjd_tdb"].to_numpy(float)
    y = sub["rv_centered"].to_numpy(float)
    dy = sub["rv_err"].to_numpy(float)

    good = np.isfinite(t) & np.isfinite(y) & np.isfinite(dy) & (dy > 0)
    t, y, dy = t[good], y[good], dy[good]

    fmin = 1.0 / P_MAX
    fmax = 1.0 / P_MIN

    ls = LombScargle(t, y, dy, fit_mean=True, center_data=True)
    freq, power = ls.autopower(
        minimum_frequency=fmin,
        maximum_frequency=fmax,
        samples_per_peak=SAMPLES_PER_PEAK,
    )

    period = 1.0 / freq
    window = classical_window(t, freq)

    order = np.argsort(period)
    return period[order], power[order], window[order]

def top_peaks(period, power, pmin, pmax, n=8):
    mask = (period >= pmin) & (period <= pmax)
    p = period[mask]
    y = power[mask]
    idx, _ = find_peaks(y)
    if len(idx) == 0:
        return []
    order = idx[np.argsort(y[idx])[::-1]]

    chosen = []
    for i in order:
        pi, yi = p[i], y[i]
        if all(abs(pi - pj) > 0.03 * pj for pj, _ in chosen):
            chosen.append((pi, yi))
        if len(chosen) >= n:
            break
    return chosen

def add_marks(ax, ymax, xmin, xmax):
    for label, p in MARKS.items():
        if p < xmin or p > xmax:
            continue
        lo, hi = p * (1 - VETO_FRAC), p * (1 + VETO_FRAC)
        ax.axvspan(lo, hi, alpha=0.10)
        ax.axvline(p, lw=0.7, alpha=0.7)
        if label in ["M10", "P53", "C70", "P_rot 36.2", "P_rot 38.7"]:
            ax.text(p, ymax * 0.98, label, rotation=90, va="top", ha="center", fontsize=7)

def panel(ax, period, series, title, xmin, xmax, ylabel, clip=True):
    ax.plot(period, series, lw=1.1)
    ax.set_xlim(xmin, xmax)

    mask = (period >= xmin) & (period <= xmax)
    y = series[mask]
    if clip and len(y) > 0:
        ymax = np.nanpercentile(y, 99.2) * 1.08
        ymin = np.nanpercentile(y, 0.5)
        ax.set_ylim(ymin, ymax)
    else:
        ymax = np.nanmax(y) if len(y) else np.nanmax(series)

    add_marks(ax, ymax, xmin, xmax)

    for p, pow_ in top_peaks(period, series, xmin, xmax, n=5):
        if pow_ <= ax.get_ylim()[1]:
            ax.plot(p, pow_, "o", ms=4)
            ax.text(p, pow_, f"{p:.2f}", fontsize=7, va="bottom", ha="left")

    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)

def main():
    df = prepare(pd.read_csv(INPUT))
    datasets = [
        ("ALL", df),
        ("ESPRESSO", df[df["instrument"].str.upper() == "ESPRESSO"].copy()),
        ("UVES", df[df["instrument"].str.upper() == "UVES"].copy()),
    ]

    ranges = [
        ("full_0p5_120", 0.5, 120.0),
        ("zoom_5_120", 5.0, 120.0),
        ("zoom_40_90", 40.0, 90.0),
    ]

    rows = []

    for tag, xmin, xmax in ranges:
        fig, axes = plt.subplots(len(datasets), 2, figsize=(17, 4.4 * len(datasets)), sharex=True)

        for r, (name, sub) in enumerate(datasets):
            period, power, window = calc(sub)

            panel(axes[r, 0], period, power, f"{name} — GLS RV", xmin, xmax, "GLS", clip=True)
            panel(axes[r, 1], period, window, f"{name} — window", xmin, xmax, "W(f)", clip=True)

            for kind, series in [("GLS", power), ("WINDOW", window)]:
                for rank, (p, val) in enumerate(top_peaks(period, series, xmin, xmax, n=10), 1):
                    rows.append({
                        "range": tag,
                        "dataset": name,
                        "kind": kind,
                        "rank": rank,
                        "period_d": p,
                        "power": val,
                    })

        for ax in axes[-1, :]:
            ax.set_xlabel("Periodo [d]")

        fig.suptitle(f"Alpha Cen B — dashboard {tag}", fontsize=15)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out = OUTDIR / f"gls_window_dashboard_{tag}.png"
        fig.savefig(out, dpi=180)
        plt.close(fig)
        print("Saved:", out)

    summary = pd.DataFrame(rows)
    summary.to_csv(OUTDIR / "dashboard_zoom_top_peaks.csv", index=False)
    print("Saved:", OUTDIR / "dashboard_zoom_top_peaks.csv")
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
