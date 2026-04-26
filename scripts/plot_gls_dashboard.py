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

VETO_PERIODS = [
    0.5,
    1.0,
    1.62,
    3.24,
    M10,
    36.2,
    38.7,
    36.2 / 2,
    38.7 / 2,
    36.2 / 3,
    38.7 / 3,
    36.2 / 4,
    38.7 / 4,
]

VETO_LABELS = [
    "0.5 d",
    "1 d",
    "1.62 d",
    "3.24 d",
    "M10",
    "P_rot 36.2",
    "P_rot 38.7",
    "P_rot/2",
    "P_rot/2",
    "P_rot/3",
    "P_rot/3",
    "P_rot/4",
    "P_rot/4",
]


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

    # Offset por instrumento
    df["rv_centered"] = df["rv"] - df.groupby("instrument")["rv"].transform("median")

    return df


def calc_gls(sub):
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

    return {
        "t": t,
        "period": period[order],
        "power": power[order],
        "window": window[order],
        "freq": freq[order],
        "ls": ls,
    }


def top_peaks(period, power, n=10, pmin=P_MIN, pmax=P_MAX):
    mask = (period >= pmin) & (period <= pmax)
    p = period[mask]
    y = power[mask]

    idx, _ = find_peaks(y)

    if len(idx) == 0:
        idx = np.arange(len(y))

    order = idx[np.argsort(y[idx])[::-1]]

    chosen = []
    for i in order:
        pi = p[i]
        yi = y[i]

        if all(abs(pi - pj) > 0.03 * pj for pj, _ in chosen):
            chosen.append((pi, yi))

        if len(chosen) >= n:
            break

    return chosen


def plot_vetos(ax, ymax):
    used = set()

    for p, lab in zip(VETO_PERIODS, VETO_LABELS):
        lo = p * (1 - VETO_FRAC)
        hi = p * (1 + VETO_FRAC)

        if hi < P_MIN or lo > P_MAX:
            continue

        ax.axvspan(lo, hi, alpha=0.12)
        ax.axvline(p, lw=0.7, alpha=0.85)

        key = round(p, 3)
        if key not in used:
            used.add(key)
            ax.text(
                p,
                ymax * 0.98,
                lab,
                rotation=90,
                va="top",
                ha="center",
                fontsize=7,
            )


def main():
    df = prepare(pd.read_csv(INPUT))

    datasets = [
        ("ALL", df),
        ("ESPRESSO", df[df["instrument"].str.upper() == "ESPRESSO"].copy()),
        ("UVES", df[df["instrument"].str.upper() == "UVES"].copy()),
    ]

    results = {}
    rows = []

    for name, sub in datasets:
        if len(sub) < 5:
            continue

        res = calc_gls(sub)
        results[name] = res

        for rank, (p, pow_) in enumerate(top_peaks(res["period"], res["power"], n=12), start=1):
            rows.append({
                "dataset": name,
                "kind": "GLS",
                "rank": rank,
                "period_d": p,
                "power": pow_,
            })

        for rank, (p, pow_) in enumerate(top_peaks(res["period"], res["window"], n=12), start=1):
            rows.append({
                "dataset": name,
                "kind": "WINDOW",
                "rank": rank,
                "period_d": p,
                "power": pow_,
            })

    summary = pd.DataFrame(rows)
    summary.to_csv(OUTDIR / "dashboard_top_peaks.csv", index=False)

    fig, axes = plt.subplots(
        nrows=len(results),
        ncols=2,
        figsize=(17, 4.7 * len(results)),
        sharex=True,
    )

    if len(results) == 1:
        axes = np.array([axes])

    for row_idx, (name, res) in enumerate(results.items()):
        period = res["period"]
        power = res["power"]
        window = res["window"]

        # GLS
        ax = axes[row_idx, 0]
        ax.plot(period, power, lw=1.2)
        ymax = np.nanmax(power)
        plot_vetos(ax, ymax)

        for p, pow_ in top_peaks(period, power, n=5):
            ax.plot(p, pow_, "o")
            ax.text(p, pow_, f"{p:.2f} d", fontsize=8, va="bottom", ha="left")

        ax.set_title(f"{name} — GLS RV")
        ax.set_ylabel("Potencia GLS")
        ax.grid(alpha=0.25)
        ax.set_xlim(P_MIN, P_MAX)

        # WINDOW
        ax = axes[row_idx, 1]
        ax.plot(period, window, lw=1.2)
        ymax = np.nanmax(window)
        plot_vetos(ax, ymax)

        for p, pow_ in top_peaks(period, window, n=5):
            ax.plot(p, pow_, "o")
            ax.text(p, pow_, f"{p:.2f} d", fontsize=8, va="bottom", ha="left")

        ax.set_title(f"{name} — classical spectral window")
        ax.set_ylabel("W(f)")
        ax.grid(alpha=0.25)
        ax.set_xlim(P_MIN, P_MAX)

    for ax in axes[-1, :]:
        ax.set_xlabel("Periodo [d]")

    fig.suptitle("Alpha Cen B — ESPRESSO vs UVES vs ALL | GLS + window corregida", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    outpng = OUTDIR / "gls_window_dashboard_all_espresso_uves.png"
    fig.savefig(outpng, dpi=180)
    plt.close(fig)

    print("Saved:")
    print(outpng)
    print(OUTDIR / "dashboard_top_peaks.csv")

    print("\n=== Top peaks summary ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
