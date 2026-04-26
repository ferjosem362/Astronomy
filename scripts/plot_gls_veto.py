#!/usr/bin/env python3
"""
plot_gls_veto.py

Periodograma GLS/Lomb-Scargle para RVs de Alpha Cen B con:
- offsets por instrumento eliminados
- bandas vetadas de armónicos de rotación
- marcas para M10, 3.24 d, 1.62 d, etc.
- spectral window function clásica:
    W(f) = |sum_j exp(2*pi*i*f*t_j)|^2 / N^2

Ejemplo:

python scripts/plot_gls_veto.py data_interim/qc/multinstrument_rvonly.csv \
  --time-col bjd_tdb \
  --rv-col rv \
  --err-col rv_err \
  --inst-col instrument \
  --outdir results/periodograms \
  --pmin 0.5 \
  --pmax 100 \
  --prot 36.2 38.7 \
  --mark 11.867564152695252 3.24 1.62
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.timeseries import LombScargle

try:
    from scipy.signal import find_peaks
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


# -----------------------------------------------------------------------------
# Column helpers
# -----------------------------------------------------------------------------
def pick_col(df, candidates, required=True):
    """
    Busca una columna por nombre exacto case-insensitive o por substring.
    """
    cols_lower = {c.lower(): c for c in df.columns}

    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]

    for c in df.columns:
        cl = c.lower()
        for cand in candidates:
            if cand.lower() in cl:
                return c

    if required:
        raise ValueError(
            f"No encuentro ninguna columna entre: {candidates}. "
            f"Columnas disponibles: {list(df.columns)}"
        )

    return None


def load_table(path):
    """
    Lee CSV/TSV/TXT intentando separador automático.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    try:
        return pd.read_csv(path, sep=None, engine="python", comment="#")
    except Exception:
        return pd.read_csv(path, delim_whitespace=True, comment="#")


# -----------------------------------------------------------------------------
# Veto / harmonics
# -----------------------------------------------------------------------------
def build_veto_periods(prot_list, extra_marks):
    """
    Construye lista de periodos vetados:
    P_rot, P_rot/2, P_rot/3, P_rot/4 + marcas extra.
    """
    veto = []
    labels = []

    for prot in prot_list:
        items = [
            (prot, f"P_rot={prot:.2f} d"),
            (prot / 2.0, f"P_rot/2={prot/2.0:.2f} d"),
            (prot / 3.0, f"P_rot/3={prot/3.0:.2f} d"),
            (prot / 4.0, f"P_rot/4={prot/4.0:.2f} d"),
        ]
        for p, lab in items:
            veto.append(p)
            labels.append(lab)

    for p in extra_marks:
        veto.append(p)
        labels.append(f"{p:.5g} d")

    return np.array(veto, dtype=float), labels


def in_veto(period, veto_centers, frac_width):
    """
    True si period cae dentro de alguna ventana relativa.
    """
    for vc in veto_centers:
        if vc <= 0:
            continue
        if abs(period - vc) <= frac_width * vc:
            return True
    return False


def nearest_veto(period, veto_centers, veto_labels):
    rows = []
    for vc, lab in zip(veto_centers, veto_labels):
        rows.append((abs(period - vc), vc, lab))
    rows.sort(key=lambda x: x[0])
    return rows[0]


# -----------------------------------------------------------------------------
# Peak finding
# -----------------------------------------------------------------------------
def unique_peak_periods(periods, power, npeaks=15, min_sep_frac=0.03):
    """
    Devuelve índices de picos principales evitando duplicados muy cercanos.
    """
    periods = np.asarray(periods)
    power = np.asarray(power)

    if HAVE_SCIPY:
        peak_idx, _ = find_peaks(power)
        if len(peak_idx) == 0:
            peak_idx = np.arange(len(power))
    else:
        peak_idx = np.arange(1, len(power) - 1)
        mask = (power[1:-1] > power[:-2]) & (power[1:-1] > power[2:])
        peak_idx = peak_idx[mask]
        if len(peak_idx) == 0:
            peak_idx = np.arange(len(power))

    order = peak_idx[np.argsort(power[peak_idx])[::-1]]

    chosen_periods = []
    chosen_idx = []

    for idx in order:
        p = periods[idx]
        if not np.isfinite(p):
            continue

        too_close = False
        for pc in chosen_periods:
            if abs(p - pc) <= min_sep_frac * pc:
                too_close = True
                break

        if not too_close:
            chosen_periods.append(p)
            chosen_idx.append(idx)

        if len(chosen_idx) >= npeaks:
            break

    return np.array(chosen_idx, dtype=int)


# -----------------------------------------------------------------------------
# Classical spectral window
# -----------------------------------------------------------------------------
def classical_window_function(t, freq):
    """
    W(f) = |sum_j exp(2*pi*i*f*t_j)|^2 / N^2

    t en días, freq en ciclos/día.
    Se resta t.min() para estabilidad numérica.
    """
    t = np.asarray(t, dtype=float)
    freq = np.asarray(freq, dtype=float)

    tt = t - np.nanmin(t)
    n = len(tt)

    # Para evitar matrices absurdamente grandes en casos enormes,
    # calculamos por chunks. Aquí N~200 y freq~miles, pero esto es robusto.
    out = np.empty_like(freq)

    chunk = 2000
    for i in range(0, len(freq), chunk):
        ff = freq[i:i + chunk]
        phase = np.exp(2j * np.pi * ff[:, None] * tt[None, :])
        out[i:i + chunk] = np.abs(np.sum(phase, axis=1)) ** 2 / (n ** 2)

    return out


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="GLS periodogram + harmonic veto + classical spectral window"
    )

    ap.add_argument("input", help="Tabla RV: csv/tsv/txt/dat")
    ap.add_argument("--outdir", default="results/periodograms")

    ap.add_argument("--time-col", default=None)
    ap.add_argument("--rv-col", default=None)
    ap.add_argument("--err-col", default=None)
    ap.add_argument("--inst-col", default=None)

    ap.add_argument("--pmin", type=float, default=0.5)
    ap.add_argument("--pmax", type=float, default=100.0)
    ap.add_argument("--samples-per-peak", type=int, default=20)

    ap.add_argument(
        "--veto-frac",
        type=float,
        default=0.03,
        help="Semiancho relativo de veto. 0.03 = ±3%%."
    )

    ap.add_argument(
        "--prot",
        nargs="*",
        type=float,
        default=[36.2, 38.7],
        help="Periodos candidatos de rotación estelar."
    )

    ap.add_argument(
        "--mark",
        nargs="*",
        type=float,
        default=[11.867564152695252, 3.24, 1.62],
        help="Periodos extra a marcar/vetar."
    )

    ap.add_argument(
        "--title",
        default="Alpha Cen B",
        help="Título base de las figuras."
    )

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_table(args.input)

    time_candidates = [
        "bjd_tdb", "bjd", "jd", "mjd", "mjd_obs",
        "time", "times", "t"
    ]
    rv_candidates = [
        "rv", "rvs", "vrad", "velocity", "vel"
    ]
    err_candidates = [
        "rv_err", "rverr", "e_rv", "erv", "sigma_rv",
        "rv_error", "err", "error"
    ]
    inst_candidates = [
        "instrument", "inst", "instr", "source", "dataset"
    ]

    tcol = args.time_col or pick_col(df, time_candidates, required=True)
    rvcol = args.rv_col or pick_col(df, rv_candidates, required=True)
    ecol = args.err_col or pick_col(df, err_candidates, required=False)
    icol = args.inst_col or pick_col(df, inst_candidates, required=False)

    needed = [tcol, rvcol]
    if ecol is not None:
        needed.append(ecol)
    if icol is not None:
        needed.append(icol)

    work = df.dropna(subset=needed).copy()

    # quality_flag = 0 si existe
    if "quality_flag" in work.columns:
        before = len(work)
        work = work[work["quality_flag"].fillna(0).astype(int) == 0].copy()
        after = len(work)
        print(f"Aplicado quality_flag==0: {before} -> {after}")

    t = np.asarray(work[tcol], dtype=float)
    y_raw = np.asarray(work[rvcol], dtype=float)

    if ecol is not None:
        dy = np.asarray(work[ecol], dtype=float)
        good_err = np.isfinite(dy) & (dy > 0)
        if not np.any(good_err):
            dy = None
        else:
            med = np.nanmedian(dy[good_err])
            dy = np.where(good_err, dy, med)
    else:
        dy = None

    # Offset por instrumento
    if icol is not None:
        work["_rv_detrended"] = (
            work[rvcol] - work.groupby(icol)[rvcol].transform("median")
        )
        y = np.asarray(work["_rv_detrended"], dtype=float)
    else:
        y = y_raw - np.nanmedian(y_raw)

    # Limpieza final
    finite = np.isfinite(t) & np.isfinite(y)
    if dy is not None:
        finite &= np.isfinite(dy) & (dy > 0)

    t = t[finite]
    y = y[finite]
    if dy is not None:
        dy = dy[finite]

    if len(t) < 5:
        raise ValueError(f"Demasiados pocos puntos válidos: N={len(t)}")

    # Frecuencias
    fmin = 1.0 / args.pmax
    fmax = 1.0 / args.pmin

    ls = LombScargle(t, y, dy, fit_mean=True, center_data=True)
    freq, power = ls.autopower(
        minimum_frequency=fmin,
        maximum_frequency=fmax,
        samples_per_peak=args.samples_per_peak,
    )
    periods = 1.0 / freq

    # Ordenar por periodo ascendente para plot legible
    order_p = np.argsort(periods)
    periods_plot = periods[order_p]
    power_plot = power[order_p]

    # Classical window usando las mismas frecuencias
    power_w = classical_window_function(t, freq)
    periods_w = 1.0 / freq

    order_w = np.argsort(periods_w)
    periods_w_plot = periods_w[order_w]
    power_w_plot = power_w[order_w]

    # FAP levels
    fap_levels = [0.1, 0.01, 0.001]
    fap_powers = []
    for fap in fap_levels:
        try:
            fap_powers.append(ls.false_alarm_level(fap))
        except Exception:
            fap_powers.append(np.nan)

    # Vetos
    veto_centers, veto_labels = build_veto_periods(args.prot, args.mark)

    # Picos GLS
    peak_idx = unique_peak_periods(periods, power, npeaks=15, min_sep_frac=0.03)

    peak_table = []
    for idx in peak_idx:
        p = float(periods[idx])
        pw = float(power[idx])
        veto_flag = in_veto(p, veto_centers, args.veto_frac)
        ndist, ncenter, nlabel = nearest_veto(p, veto_centers, veto_labels)
        peak_table.append({
            "period_d": p,
            "power": pw,
            "veto": veto_flag,
            "nearest_veto": nlabel,
            "nearest_veto_period": ncenter,
            "delta_to_nearest_veto_d": ndist,
        })

    peak_table = sorted(peak_table, key=lambda r: r["power"], reverse=True)

    # Picos window
    wpeak_idx = unique_peak_periods(periods_w, power_w, npeaks=15, min_sep_frac=0.03)
    wpeak_table = []
    for idx in wpeak_idx:
        p = float(periods_w[idx])
        pw = float(power_w[idx])
        veto_flag = in_veto(p, veto_centers, args.veto_frac)
        ndist, ncenter, nlabel = nearest_veto(p, veto_centers, veto_labels)
        wpeak_table.append({
            "period_d": p,
            "window_power": pw,
            "veto": veto_flag,
            "nearest_veto": nlabel,
            "nearest_veto_period": ncenter,
            "delta_to_nearest_veto_d": ndist,
        })

    wpeak_table = sorted(wpeak_table, key=lambda r: r["window_power"], reverse=True)

    # -------------------------------------------------------------------------
    # FIG 1: GLS
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 7))

    ax.plot(periods_plot, power_plot, lw=1.4)
    ax.set_xlim(args.pmin, args.pmax)
    ax.set_xlabel("Periodo [d]")
    ax.set_ylabel("Potencia GLS")
    ax.set_title(f"{args.title} — periodograma GLS con bandas vetadas")
    ax.grid(alpha=0.25)

    ymax = np.nanmax(power_plot)
    ymin = np.nanmin(power_plot)

    used = set()
    for vc, lab in zip(veto_centers, veto_labels):
        if not np.isfinite(vc) or vc <= 0:
            continue

        lo = vc * (1.0 - args.veto_frac)
        hi = vc * (1.0 + args.veto_frac)

        if hi < args.pmin or lo > args.pmax:
            continue

        ax.axvspan(lo, hi, alpha=0.13)
        ax.axvline(vc, lw=0.8, alpha=0.9)

        key = round(vc, 4)
        if key not in used:
            used.add(key)
            ax.text(
                vc,
                ymax * 0.98,
                lab,
                rotation=90,
                va="top",
                ha="center",
                fontsize=8,
            )

    # M10 destacado
    m10 = 11.867564152695252
    if args.pmin <= m10 <= args.pmax:
        ax.axvline(m10, lw=2.0)
        ax.text(
            m10,
            ymax * 0.74,
            "M10 = 11.8676 d",
            rotation=90,
            va="top",
            ha="left",
            fontsize=9,
        )

    for lev, fp in zip(fap_levels, fap_powers):
        if np.isfinite(fp):
            ax.axhline(fp, ls="--", lw=1.0)
            ax.text(
                args.pmax * 0.98,
                fp,
                f"FAP={lev:g}",
                va="bottom",
                ha="right",
                fontsize=8,
            )

    # Etiquetar top peaks
    for row in peak_table[:10]:
        p = row["period_d"]
        pw = row["power"]
        veto_flag = row["veto"]

        if not (args.pmin <= p <= args.pmax):
            continue

        ax.plot([p], [pw], marker="o")
        txt = f"{p:.2f} d"
        if veto_flag:
            txt += " [VETO]"

        ax.text(
            p,
            pw,
            txt,
            fontsize=8,
            va="bottom",
            ha="left",
        )

    f_gls = outdir / "gls_periodogram_veto.png"
    fig.tight_layout()
    fig.savefig(f_gls, dpi=180)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # FIG 2: Classical window
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 6))

    ax.plot(periods_w_plot, power_w_plot, lw=1.4)
    ax.set_xlim(args.pmin, args.pmax)
    ax.set_xlabel("Periodo [d]")
    ax.set_ylabel("W(f)")
    ax.set_title(f"{args.title} — classical spectral window")
    ax.grid(alpha=0.25)

    ymax_w = np.nanmax(power_w_plot)

    used = set()
    for vc, lab in zip(veto_centers, veto_labels):
        if not np.isfinite(vc) or vc <= 0:
            continue

        lo = vc * (1.0 - args.veto_frac)
        hi = vc * (1.0 + args.veto_frac)

        if hi < args.pmin or lo > args.pmax:
            continue

        ax.axvspan(lo, hi, alpha=0.13)
        ax.axvline(vc, lw=0.8, alpha=0.9)

        key = round(vc, 4)
        if key not in used:
            used.add(key)
            ax.text(
                vc,
                ymax_w * 0.98,
                lab,
                rotation=90,
                va="top",
                ha="center",
                fontsize=8,
            )

    for row in wpeak_table[:10]:
        p = row["period_d"]
        pw = row["window_power"]

        if not (args.pmin <= p <= args.pmax):
            continue

        ax.plot([p], [pw], marker="o")
        ax.text(
            p,
            pw,
            f"{p:.2f} d",
            fontsize=8,
            va="bottom",
            ha="left",
        )

    f_window = outdir / "window_function.png"
    fig.tight_layout()
    fig.savefig(f_window, dpi=180)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # FIG 3: GLS + window overlay normalizado
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(13, 6))

    pwr_norm = power_plot / np.nanmax(power_plot)
    win_norm = power_w_plot / np.nanmax(power_w_plot)

    ax.plot(periods_plot, pwr_norm, lw=1.3, label="GLS RV normalizado")
    ax.plot(periods_w_plot, win_norm, lw=1.2, alpha=0.75, label="Window normalizada")

    ax.set_xlim(args.pmin, args.pmax)
    ax.set_xlabel("Periodo [d]")
    ax.set_ylabel("Potencia normalizada")
    ax.set_title(f"{args.title} — GLS vs window function")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    for vc, lab in zip(veto_centers, veto_labels):
        if not np.isfinite(vc) or vc <= 0:
            continue

        lo = vc * (1.0 - args.veto_frac)
        hi = vc * (1.0 + args.veto_frac)

        if hi < args.pmin or lo > args.pmax:
            continue

        ax.axvspan(lo, hi, alpha=0.10)
        ax.axvline(vc, lw=0.7, alpha=0.8)

    f_overlay = outdir / "gls_vs_window_overlay.png"
    fig.tight_layout()
    fig.savefig(f_overlay, dpi=180)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Save tables
    # -------------------------------------------------------------------------
    peaks_csv = outdir / "gls_top_peaks.csv"
    window_csv = outdir / "window_top_peaks.csv"

    pd.DataFrame(peak_table).to_csv(peaks_csv, index=False)
    pd.DataFrame(wpeak_table).to_csv(window_csv, index=False)

    # -------------------------------------------------------------------------
    # Terminal report
    # -------------------------------------------------------------------------
    print("\n=== Input ===")
    print(args.input)
    print(f"N used: {len(t)}")

    print("\n=== Columnas usadas ===")
    print("time :", tcol)
    print("rv   :", rvcol)
    print("err  :", ecol)
    print("inst :", icol)

    if icol is not None:
        print("\n=== Puntos por instrumento ===")
        print(work.loc[finite, icol].value_counts().to_string())

    print("\n=== Top GLS peaks ===")
    for i, row in enumerate(peak_table[:15], start=1):
        print(
            f"{i:2d}. "
            f"P={row['period_d']:10.5f} d   "
            f"power={row['power']:10.6f}   "
            f"veto={str(row['veto']):5s}   "
            f"nearest={row['nearest_veto']} "
            f"(Δ={row['delta_to_nearest_veto_d']:.5f} d)"
        )

    print("\n=== Top WINDOW peaks ===")
    for i, row in enumerate(wpeak_table[:15], start=1):
        print(
            f"{i:2d}. "
            f"P={row['period_d']:10.5f} d   "
            f"W={row['window_power']:10.6f}   "
            f"veto={str(row['veto']):5s}   "
            f"nearest={row['nearest_veto']} "
            f"(Δ={row['delta_to_nearest_veto_d']:.5f} d)"
        )

    print("\n=== Distancias de M10 a armónicos/marcas vetadas ===")
    rows = []
    for vc, lab in zip(veto_centers, veto_labels):
        dist = abs(m10 - vc)
        frac = dist / vc if vc != 0 else np.nan
        rows.append((lab, vc, dist, frac))

    rows = sorted(rows, key=lambda x: x[2])
    for lab, vc, dist, frac in rows:
        print(
            f"{lab:18s}  "
            f"center={vc:10.5f} d  "
            f"|ΔP|={dist:8.5f} d  "
            f"frac={100*frac:7.3f}%"
        )

    print("\n=== Outputs ===")
    print(f_gls)
    print(f_window)
    print(f_overlay)
    print(peaks_csv)
    print(window_csv)


if __name__ == "__main__":
    main()
