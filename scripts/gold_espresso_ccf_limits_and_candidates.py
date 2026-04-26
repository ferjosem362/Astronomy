#!/usr/bin/env python3
"""
Gold ESPRESSO CCF planet search + detection limits.

Inputs expected:
  data_interim/qc/espresso_eso_ccf_curated_hr.csv

Outputs:
  results/gold_espresso_ccf/
    00_clean_input.csv
    01_segment_detrended.csv
    02_gls_raw.csv/.png
    03_gls_detrended.csv/.png
    04_window.csv/.png
    05_activity_periodograms.csv/.png
    06_candidate_alias_activity_audit.csv
    07_injection_recovery_matrix.csv
    08_detection_limits_mass.csv
    09_detection_limits_map.png
    10_detection_limits_slices.png
    11_summary.txt

Parallel injection recovery:
  Use --job-index and --n-jobs, then merge with --merge-only.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import warnings
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import signal
from scipy.optimize import curve_fit
from astropy.timeseries import LombScargle

warnings.filterwarnings("ignore", category=RuntimeWarning)

DAY_ALIASES = np.array([0.5, 0.970108077359177, 0.9941836673034685, 1.0, 1.008383808261102, 1.028872682954664])
TARGET_PERIODS = np.array([0.5, 0.970108077359177, 1.0, 1.028872682954664, 2.024, 2.257, 3.236,
                           5.0, 7.0, 16.3009, 18.9805, 26.6578, 32.3315, 35.6725,
                           46.8600, 55.2678, 89.4351, 123.3857, 170.4414, 365.25])


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def first_col(df: pd.DataFrame, candidates: list[str], required=True) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(f"No encuentro ninguna columna de {candidates}. Columnas: {list(df.columns)}")
    return None


def robust_mad(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan
    med = np.median(x)
    return 1.4826 * np.median(np.abs(x - med))


def sigma_clip_df(df, col, nsig=6.0, group_col=None):
    if group_col is None:
        med = np.nanmedian(df[col])
        mad = robust_mad(df[col])
        if not np.isfinite(mad) or mad <= 0:
            return df.copy()
        return df[np.abs(df[col] - med) <= nsig * mad].copy()

    keep = []
    for _, g in df.groupby(group_col):
        med = np.nanmedian(g[col])
        mad = robust_mad(g[col])
        if not np.isfinite(mad) or mad <= 0:
            keep.append(g)
        else:
            keep.append(g[np.abs(g[col] - med) <= nsig * mad])
    return pd.concat(keep, ignore_index=True)


def load_clean(input_path: Path) -> pd.DataFrame:
    df = pd.read_csv(input_path)

    time_col = first_col(df, ["time", "bjd", "bjd_approx", "bjd_dace"])
    rv_col = first_col(df, ["rv", "rv_mps", "rv_centered", "rv_centered_by_model_inst"])
    err_col = first_col(df, ["rv_err", "rv_err_mps", "e_rv", "err"])
    inst_col = first_col(df, ["model_inst", "sub_inst", "inst", "instrument"], required=False)
    mjd_col = first_col(df, ["mjd_obs", "mjd", "night_mjd"], required=False)

    out = pd.DataFrame()
    out["time"] = pd.to_numeric(df[time_col], errors="coerce")
    out["rv"] = pd.to_numeric(df[rv_col], errors="coerce")
    out["rv_err"] = pd.to_numeric(df[err_col], errors="coerce")

    if mjd_col:
        out["mjd_obs"] = pd.to_numeric(df[mjd_col], errors="coerce")
    else:
        out["mjd_obs"] = out["time"] - 2400000.5

    if inst_col:
        out["model_inst"] = df[inst_col].astype(str).fillna("ESPRESSO")
    else:
        out["model_inst"] = "ESPRESSO"

    for src, dst in [
        ("berv", "berv"), ("berv_mps", "berv"),
        ("ccf_fwhm", "ccf_fwhm"), ("ccf_fwhm_mps", "ccf_fwhm"),
        ("ccf_bispan", "ccf_bispan"), ("ccf_bis_span_mps", "ccf_bispan"),
        ("ccf_contrast", "ccf_contrast"),
        ("mask", "mask"), ("ccf_mask", "mask"),
        ("ins_mode", "ins_mode"),
        ("filename", "filename"),
    ]:
        if src in df.columns and dst not in out.columns:
            out[dst] = df[src]

    numeric_aux = ["berv", "ccf_fwhm", "ccf_bispan", "ccf_contrast"]
    for c in numeric_aux:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")

    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["time", "rv", "rv_err"])
    out = out[out["rv_err"] > 0].copy()
    out = out.sort_values("time").reset_index(drop=True)

    # safety: if errors are clearly in mm/s, convert to m/s. Existing curated should already be m/s.
    mederr = np.nanmedian(out["rv_err"])
    if mederr > 10:
        out["rv_err"] = out["rv_err"] / 1000.0
        out["rv_err_unit_fixed"] = True
    else:
        out["rv_err_unit_fixed"] = False

    out["night_mjd"] = np.floor(out["mjd_obs"]).astype(int)

    # remove obvious non-HR/odd mask rows if present
    if "ins_mode" in out.columns:
        out = out[out["ins_mode"].astype(str).str.contains("HR", case=False, na=False)].copy()

    if "mask" in out.columns:
        # keep K2 and unknown; drop G2 if mixed
        mask_s = out["mask"].astype(str)
        if (mask_s == "K2").sum() > 0:
            out = out[(mask_s == "K2") | mask_s.isin(["nan", "None", ""])].copy()

    # remove pathological RV chunks: alpha Cen B HR should be around -26 km/s after CCF selection
    # keep this conservative cut to remove known ~-18.9 km/s bad/SOBF contamination blocks.
    out = out[(out["rv"] < -24000) & (out["rv"] > -28000)].copy()

    out = sigma_clip_df(out, "rv", nsig=10.0, group_col="model_inst")
    out = out.reset_index(drop=True)
    return out


def weighted_poly_detrend_segment(t, y, e, degree):
    t0 = np.nanmedian(t)
    x = t - t0
    w = 1.0 / np.maximum(e, 1e-6)
    deg = min(degree, max(0, len(t) - 2))
    coef = np.polyfit(x, y, deg=deg, w=w)
    trend = np.polyval(coef, x)
    return trend, coef, t0


def detrend_by_segment(df: pd.DataFrame, degree_long=2, degree_short=1) -> tuple[pd.DataFrame, pd.DataFrame]:
    chunks = []
    rows = []
    for inst, g in df.groupby("model_inst", sort=False):
        g = g.sort_values("time").copy()
        baseline = g["time"].max() - g["time"].min()
        deg = degree_long if baseline > 100 else degree_short
        trend, coef, t0 = weighted_poly_detrend_segment(
            g["time"].to_numpy(float),
            g["rv"].to_numpy(float),
            g["rv_err"].to_numpy(float),
            deg,
        )
        g["trend_poly"] = trend
        g["rv_segment_centered"] = g["rv"] - np.average(g["rv"], weights=1 / g["rv_err"]**2)
        g["rv_detrended"] = g["rv"] - trend
        g["rv_detrended"] -= np.average(g["rv_detrended"], weights=1 / g["rv_err"]**2)
        chunks.append(g)
        rows.append({
            "model_inst": inst,
            "n": len(g),
            "time_min": g["time"].min(),
            "time_max": g["time"].max(),
            "baseline_d": baseline,
            "degree": deg,
            "rv_std": np.nanstd(g["rv"], ddof=1),
            "rv_segment_centered_std": np.nanstd(g["rv_segment_centered"], ddof=1),
            "rv_detrended_std": np.nanstd(g["rv_detrended"], ddof=1),
            "rv_err_median": np.nanmedian(g["rv_err"]),
            "poly_t0": t0,
            "poly_coef_json": json.dumps([float(x) for x in np.atleast_1d(coef)]),
        })
    out = pd.concat(chunks, ignore_index=True).sort_values("time").reset_index(drop=True)
    summ = pd.DataFrame(rows)
    return out, summ


def make_frequency_grid(t, min_period=0.5, max_period=None, samples_per_peak=20):
    baseline = float(np.max(t) - np.min(t))
    if max_period is None:
        max_period = 0.9 * baseline
    min_freq = 1.0 / max_period
    max_freq = 1.0 / min_period
    dfreq = 1.0 / (baseline * samples_per_peak)
    n = int(np.ceil((max_freq - min_freq) / dfreq))
    freq = min_freq + dfreq * np.arange(n + 1)
    freq = freq[freq > 0]
    return freq


def gls(t, y, e=None, min_period=0.5, max_period=None, samples_per_peak=20):
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    if e is None:
        e = np.ones_like(y)
    else:
        e = np.asarray(e, float)
    freq = make_frequency_grid(t, min_period=min_period, max_period=max_period, samples_per_peak=samples_per_peak)
    ls = LombScargle(t, y, e, fit_mean=True, center_data=True)
    power = ls.power(freq, method="auto")
    period = 1.0 / freq
    order = np.argsort(period)
    res = pd.DataFrame({"period_d": period[order], "frequency_1_d": freq[order], "power": power[order]})
    return res, ls


def top_independent_peaks(pg: pd.DataFrame, n=30, exclusion_frac=0.03):
    d = pg.sort_values("power", ascending=False).copy()
    chosen = []
    for _, row in d.iterrows():
        p = float(row["period_d"])
        if all(abs(p - c["period_d"]) / p > exclusion_frac for c in chosen):
            chosen.append(row.to_dict())
        if len(chosen) >= n:
            break
    return pd.DataFrame(chosen)


def plot_periodogram(pg, top, out_png, title, target_periods=TARGET_PERIODS):
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.plot(pg["period_d"], pg["power"], lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Period [d]")
    ax.set_ylabel("GLS power")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    for p in target_periods:
        if pg["period_d"].min() <= p <= pg["period_d"].max():
            ax.axvline(p, ls="--", lw=0.7, alpha=0.35)
    if top is not None and len(top):
        for _, r in top.head(10).iterrows():
            ax.annotate(f"{r['period_d']:.2f} d", xy=(r["period_d"], r["power"]),
                        xytext=(5, 5), textcoords="offset points", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def plot_timeseries(df, out_png):
    fig, ax = plt.subplots(figsize=(14, 6))
    for inst, g in df.groupby("model_inst"):
        ax.errorbar(g["time"], g["rv_detrended"], yerr=g["rv_err"], fmt=".", ms=3, alpha=0.75, label=inst)
    ax.set_xlabel("BJD")
    ax.set_ylabel("Detrended RV [m/s]")
    ax.set_title("ESPRESSO CCF gold detrended RV")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def window_periodogram(t, min_period=0.5, max_period=None, samples_per_peak=20):
    y = np.ones_like(np.asarray(t, float))
    e = np.ones_like(y)
    freq = make_frequency_grid(t, min_period=min_period, max_period=max_period, samples_per_peak=samples_per_peak)
    # Centering all-ones gives zero signal, so use spectral window explicitly.
    tau = t - np.min(t)
    z = np.exp(2j * np.pi * np.outer(freq, tau)).sum(axis=1)
    power = (np.abs(z) / len(t)) ** 2
    period = 1.0 / freq
    order = np.argsort(period)
    return pd.DataFrame({"period_d": period[order], "frequency_1_d": freq[order], "window_power": power[order]})


def local_power(pg, p, col="power", frac=0.03):
    d = pg[(pg["period_d"] >= p * (1 - frac)) & (pg["period_d"] <= p * (1 + frac))]
    if len(d) == 0:
        idx = (pg["period_d"] - p).abs().idxmin()
        return float(pg.loc[idx, "period_d"]), float(pg.loc[idx, col])
    row = d.loc[d[col].idxmax()]
    return float(row["period_d"]), float(row[col])


def audit_candidates(pg_det, pg_win, df, top, out_csv):
    rows = []
    activity_cols = [c for c in ["berv", "ccf_fwhm", "ccf_bispan", "ccf_contrast"] if c in df.columns and np.isfinite(pd.to_numeric(df[c], errors="coerce")).sum() > 20]

    activity_pgs = {}
    for c in activity_cols:
        y = pd.to_numeric(df[c], errors="coerce").to_numpy(float)
        ok = np.isfinite(y)
        if ok.sum() > 20:
            activity_pgs[c], _ = gls(df.loc[ok, "time"], y[ok], None,
                                     min_period=pg_det["period_d"].min(),
                                     max_period=pg_det["period_d"].max(),
                                     samples_per_peak=10)

    for _, r in top.iterrows():
        p = float(r["period_d"])
        _, wloc = local_power(pg_win.rename(columns={"window_power": "power"}), p, col="power", frac=0.02)

        day_alias_dist = np.min(np.abs(p - DAY_ALIASES) / p)
        alias_flag = day_alias_dist < 0.03

        row = {
            "period_d": p,
            "rv_power": float(r["power"]),
            "window_power_near": wloc,
            "rv_over_window": float(r["power"]) / max(wloc, 1e-12),
            "near_daily_alias": bool(alias_flag),
            "nearest_daily_alias_d": float(DAY_ALIASES[np.argmin(np.abs(p - DAY_ALIASES))]),
        }

        for c, apg in activity_pgs.items():
            _, apow = local_power(apg, p, col="power", frac=0.03)
            row[f"{c}_power_near"] = apow
            row[f"rv_over_{c}"] = float(r["power"]) / max(apow, 1e-12)

        rows.append(row)

    aud = pd.DataFrame(rows)
    aud.to_csv(out_csv, index=False)
    return aud, activity_pgs


def plot_activity_periodograms(activity_pgs, pg_det, out_png):
    if not activity_pgs:
        return
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.plot(pg_det["period_d"], pg_det["power"] / np.nanmax(pg_det["power"]), lw=1.2, label="RV detrended")
    for c, pg in activity_pgs.items():
        ax.plot(pg["period_d"], pg["power"] / np.nanmax(pg["power"]), lw=0.9, alpha=0.8, label=c)
    ax.set_xscale("log")
    ax.set_xlabel("Period [d]")
    ax.set_ylabel("Normalized GLS power")
    ax.set_title("RV vs activity/diagnostic periodograms")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def sinusoid_design(t, periods, include_const=True):
    cols = []
    names = []
    if include_const:
        cols.append(np.ones_like(t))
        names.append("const")
    t0 = np.median(t)
    x = t - t0
    for p in periods:
        w = 2 * np.pi / p
        cols.append(np.sin(w * x))
        names.append(f"sin_{p:.8g}")
        cols.append(np.cos(w * x))
        names.append(f"cos_{p:.8g}")
    return np.vstack(cols).T, names


def weighted_linear_fit(t, y, e, periods):
    X, names = sinusoid_design(t, periods, include_const=True)
    w = 1.0 / np.maximum(e, 1e-6)
    Xw = X * w[:, None]
    yw = y * w
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    model = X @ beta
    resid = y - model
    chi2 = np.sum((resid / e) ** 2)
    k = len(beta)
    n = len(y)
    bic = chi2 + k * np.log(n)
    rms = np.std(resid, ddof=1)
    wrms = np.sqrt(np.average(resid**2, weights=1/e**2))
    amps = {}
    for p in periods:
        sname = f"sin_{p:.8g}"
        cname = f"cos_{p:.8g}"
        s = beta[names.index(sname)]
        c = beta[names.index(cname)]
        amps[f"amp_{p:.8g}"] = float(np.hypot(s, c))
        amps[f"phase_{p:.8g}"] = float(np.arctan2(c, s))
    return {"chi2": chi2, "bic": bic, "rms": rms, "wrms": wrms, "k": k, "n": n, **amps}, model, resid


def candidate_model_grid(df, candidate_periods, out_csv):
    t = df["time"].to_numpy(float)
    y = df["rv_detrended"].to_numpy(float)
    e = df["rv_err"].to_numpy(float)

    period_sets = [
        [],
        [35.6724877],
        [32.3314560],
        [46.8600303],
        [55.2678282],
        [89.4351158],
        [123.385733],
        [170.441409],
        [35.6724877, 32.3314560],
        [35.6724877, 46.8600303],
        [35.6724877, 32.3314560, 46.8600303],
        [0.970108077, 35.6724877],
        [1.028872683, 35.6724877],
        [16.3009094, 18.9804574, 35.6724877],
    ]

    # Add top periods automatically, excluding daily aliases.
    for p in candidate_periods[:8]:
        if np.min(np.abs(p - DAY_ALIASES) / p) > 0.03:
            period_sets.append([float(p)])

    seen = set()
    rows = []
    models = {}
    for ps in period_sets:
        key = tuple(np.round(ps, 8))
        if key in seen:
            continue
        seen.add(key)
        fit, model, resid = weighted_linear_fit(t, y, e, ps)
        row = {
            "model": "NULL_offsets_only" if len(ps) == 0 else "P_" + "_".join(f"{p:.4g}" for p in ps),
            "periods": ",".join(map(str, ps)),
            "n_periods": len(ps),
            **fit,
        }
        rows.append(row)
        models[row["model"]] = (model, resid)

    res = pd.DataFrame(rows).sort_values("bic").reset_index(drop=True)
    best_bic = res["bic"].min()
    null_bic = float(res.loc[res["n_periods"].eq(0), "bic"].iloc[0])
    res["delta_bic_vs_best"] = res["bic"] - best_bic
    res["delta_bic_vs_null"] = res["bic"] - null_bic
    res.to_csv(out_csv, index=False)

    best_name = res.iloc[0]["model"]
    best_model, best_resid = models[best_name]
    out = df.copy()
    out["best_model_rv"] = best_model
    out["best_model_resid"] = best_resid
    return res, out


def mass_to_k_mps(m_earth, period_d, mstar_msun=0.934):
    # Circular orbit, Mp << Mstar, sin i = 1.
    p_yr = period_d / 365.25
    return 0.08945 * m_earth * (p_yr ** (-1/3)) * (mstar_msun ** (-2/3))


def k_to_mass_earth(k_mps, period_d, mstar_msun=0.934):
    p_yr = period_d / 365.25
    return k_mps / (0.08945 * (p_yr ** (-1/3)) * (mstar_msun ** (-2/3)))


def inject_signal(t, period, k, phase):
    return k * np.sin(2 * np.pi * t / period + phase)


def recovery_test_one(t, y, e, period, mass, phase, baseline, power_threshold, mstar):
    k = mass_to_k_mps(mass, period, mstar)
    yinj = y + inject_signal(t, period, k, phase)
    pg, _ = gls(t, yinj, e, min_period=0.5, max_period=0.9 * baseline, samples_per_peak=12)

    local = pg[(pg["period_d"] >= period * 0.98) & (pg["period_d"] <= period * 1.02)]
    if len(local) == 0:
        idx = (pg["period_d"] - period).abs().idxmin()
        best_p = float(pg.loc[idx, "period_d"])
        best_power = float(pg.loc[idx, "power"])
    else:
        row = local.loc[local["power"].idxmax()]
        best_p = float(row["period_d"])
        best_power = float(row["power"])

    recovered_period = abs(best_p - period) / period <= 0.02
    recovered_power = best_power >= power_threshold
    return recovered_period and recovered_power, best_p, best_power, k


def run_injection_grid(df, outdir, periods, masses, n_phase, job_index, n_jobs, seed, mstar, power_threshold_quantile=0.999):
    rng = np.random.default_rng(seed + job_index)
    t = df["time"].to_numpy(float)
    y = df["rv_detrended"].to_numpy(float)
    e = df["rv_err"].to_numpy(float)
    baseline = float(np.max(t) - np.min(t))

    null_powers = []
    for _ in range(200):
        perm = rng.permutation(y)
        pg, _ = gls(t, perm, e, min_period=0.5, max_period=0.9 * baseline, samples_per_peak=8)
        null_powers.append(pg["power"].max())
    power_threshold = float(np.quantile(null_powers, power_threshold_quantile))

    cells = []
    for ip, p in enumerate(periods):
        for im, m in enumerate(masses):
            cells.append((ip, im, float(p), float(m)))

    cells = [c for i, c in enumerate(cells) if i % n_jobs == job_index]

    rows = []
    total = len(cells)
    for j, (ip, im, p, m) in enumerate(cells, start=1):
        successes = 0
        best_periods = []
        best_powers = []
        ks = []
        phases = rng.uniform(0, 2*np.pi, size=n_phase)
        for ph in phases:
            ok, bp, bpow, k = recovery_test_one(t, y, e, p, m, ph, baseline, power_threshold, mstar)
            successes += int(ok)
            best_periods.append(bp)
            best_powers.append(bpow)
            ks.append(k)
        rows.append({
            "period_d": p,
            "mass_earth": m,
            "k_mps": float(np.median(ks)),
            "n_trials": n_phase,
            "n_recovered": successes,
            "recovery_rate": successes / n_phase,
            "median_best_period_d": float(np.median(best_periods)),
            "median_best_power": float(np.median(best_powers)),
            "power_threshold": power_threshold,
            "job_index": job_index,
            "n_jobs": n_jobs,
        })
        if j % 10 == 0 or j == total:
            print(f"[job {job_index}] cell {j}/{total}: P={p:.3f} M={m:.3f} rec={successes}/{n_phase}", flush=True)

    part = pd.DataFrame(rows)
    part_path = outdir / f"07_injection_recovery_matrix.part{job_index:04d}.csv"
    part.to_csv(part_path, index=False)
    print("Saved:", part_path)
    return part_path


def merge_injection_parts(outdir):
    parts = sorted(outdir.glob("07_injection_recovery_matrix.part*.csv"))
    if not parts:
        raise FileNotFoundError("No hay parts de injection recovery para mergear.")
    df = pd.concat([pd.read_csv(p) for p in parts], ignore_index=True)
    df = df.sort_values(["period_d", "mass_earth"]).reset_index(drop=True)
    out = outdir / "07_injection_recovery_matrix.csv"
    df.to_csv(out, index=False)
    print("Saved:", out)
    return df


def detection_limits_from_matrix(mat, outdir):
    rows = []
    for p, g in mat.groupby("period_d"):
        g = g.sort_values("mass_earth")
        hit = g[g["recovery_rate"] >= 0.99]
        if len(hit):
            lim = hit.iloc[0]
            rows.append({
                "period_d": p,
                "mass_earth_99": lim["mass_earth"],
                "k_mps_99": lim["k_mps"],
                "recovery_at_limit": lim["recovery_rate"],
            })
        else:
            rows.append({
                "period_d": p,
                "mass_earth_99": np.nan,
                "k_mps_99": np.nan,
                "recovery_at_limit": np.nan,
            })
    lims = pd.DataFrame(rows)
    out = outdir / "08_detection_limits_mass.csv"
    lims.to_csv(out, index=False)
    print("Saved:", out)
    return lims


def plot_detection_maps(mat, lims, outdir):
    pivot = mat.pivot_table(index="mass_earth", columns="period_d", values="recovery_rate", aggfunc="mean")
    periods = pivot.columns.to_numpy(float)
    masses = pivot.index.to_numpy(float)
    Z = pivot.to_numpy(float) * 100

    fig, ax = plt.subplots(figsize=(14, 7))
    im = ax.pcolormesh(periods, masses, Z, shading="auto", vmin=0, vmax=100)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Period [d]")
    ax.set_ylabel("M sin i [Mearth]")
    ax.set_title("ESPRESSO CCF gold injection-recovery detection limits")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Recovery rate [%]")

    if lims is not None and len(lims):
        ok = lims["mass_earth_99"].notna()
        ax.plot(lims.loc[ok, "period_d"], lims.loc[ok, "mass_earth_99"], lw=2.5, label="99% recovery")
        ax.legend()

    fig.tight_layout()
    fig.savefig(outdir / "09_detection_limits_map.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 6))
    selected = [3, 5, 10, 30, 53, 100, 200, 500, 1000]
    for p0 in selected:
        p = periods[np.argmin(np.abs(periods - p0))]
        g = mat[np.isclose(mat["period_d"], p)].sort_values("mass_earth")
        ax.plot(g["mass_earth"], 100*g["recovery_rate"], marker="o", label=f"P={p:.1f} d")
    ax.axhline(99, ls="--", lw=1)
    ax.axhline(50, ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("M sin i [Mearth]")
    ax.set_ylabel("Recovery rate [%]")
    ax.set_title("Injection-recovery slices")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "10_detection_limits_slices.png", dpi=200)
    plt.close(fig)


def write_summary(outdir, clean, det_summ, top_raw, top_det, audit, model_summary, lims=None):
    p = outdir / "11_summary.txt"
    with open(p, "w") as f:
        f.write("GOLD ESPRESSO CCF SUMMARY\n")
        f.write("="*80 + "\n\n")
        f.write(f"N clean rows: {len(clean)}\n")
        f.write(f"Time baseline d: {clean['time'].max() - clean['time'].min():.3f}\n")
        f.write(f"RV err median m/s: {clean['rv_err'].median():.4f}\n\n")

        f.write("DETREND SUMMARY\n")
        f.write(det_summ.to_string(index=False))
        f.write("\n\nTOP RAW GLS\n")
        f.write(top_raw.head(20).to_string(index=False))
        f.write("\n\nTOP DETRENDED GLS\n")
        f.write(top_det.head(30).to_string(index=False))
        f.write("\n\nCANDIDATE AUDIT\n")
        f.write(audit.head(30).to_string(index=False))
        f.write("\n\nMODEL SUMMARY\n")
        f.write(model_summary.head(30).to_string(index=False))
        if lims is not None:
            f.write("\n\nDETECTION LIMITS\n")
            f.write(lims.to_string(index=False))
        f.write("\n\nDONE\n")
    print("Saved:", p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="data_interim/qc/espresso_eso_ccf_curated_hr.csv")
    ap.add_argument("--outdir", default="results/gold_espresso_ccf")
    ap.add_argument("--skip-injections", action="store_true")
    ap.add_argument("--merge-only", action="store_true")
    ap.add_argument("--job-index", type=int, default=0)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--n-phase", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260425)
    ap.add_argument("--mstar", type=float, default=0.934)
    ap.add_argument("--period-grid", default="")
    ap.add_argument("--mass-grid", default="")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    ensure_dir(outdir)

    if args.merge_only:
        mat = merge_injection_parts(outdir)
        lims = detection_limits_from_matrix(mat, outdir)
        plot_detection_maps(mat, lims, outdir)
        return

    input_path = Path(args.input)
    clean = load_clean(input_path)
    clean.to_csv(outdir / "00_clean_input.csv", index=False)
    print("Input:", input_path)
    print("Clean shape:", clean.shape)
    print("Saved:", outdir / "00_clean_input.csv")

    det, det_summ = detrend_by_segment(clean)
    det.to_csv(outdir / "01_segment_detrended.csv", index=False)
    det_summ.to_csv(outdir / "01_segment_detrended_summary.csv", index=False)
    print("Saved:", outdir / "01_segment_detrended.csv")
    print("Saved:", outdir / "01_segment_detrended_summary.csv")

    plot_timeseries(det, outdir / "01_timeseries_detrended.png")

    baseline = det["time"].max() - det["time"].min()
    max_period = 0.9 * baseline

    pg_raw, _ = gls(det["time"], det["rv_segment_centered"], det["rv_err"], min_period=0.5, max_period=max_period, samples_per_peak=25)
    pg_raw.to_csv(outdir / "02_gls_raw_segment_centered.csv", index=False)
    top_raw = top_independent_peaks(pg_raw, n=40)
    top_raw.to_csv(outdir / "02_gls_raw_segment_centered_top_peaks.csv", index=False)
    plot_periodogram(pg_raw, top_raw, outdir / "02_gls_raw_segment_centered.png", "Raw segment-centered GLS")

    pg_det, _ = gls(det["time"], det["rv_detrended"], det["rv_err"], min_period=0.5, max_period=max_period, samples_per_peak=25)
    pg_det.to_csv(outdir / "03_gls_detrended.csv", index=False)
    top_det = top_independent_peaks(pg_det, n=50)
    top_det.to_csv(outdir / "03_gls_detrended_top_peaks.csv", index=False)
    plot_periodogram(pg_det, top_det, outdir / "03_gls_detrended.png", "Detrended GLS")

    pg_win = window_periodogram(det["time"].to_numpy(float), min_period=0.5, max_period=max_period, samples_per_peak=25)
    pg_win.to_csv(outdir / "04_window.csv", index=False)
    tmp = pg_win.rename(columns={"window_power": "power"})
    top_win = top_independent_peaks(tmp, n=40)
    top_win.to_csv(outdir / "04_window_top_peaks.csv", index=False)
    plot_periodogram(tmp, top_win, outdir / "04_window.png", "Spectral window")

    audit, activity_pgs = audit_candidates(pg_det, pg_win, det, top_det, outdir / "06_candidate_alias_activity_audit.csv")
    plot_activity_periodograms(activity_pgs, pg_det, outdir / "05_activity_periodograms.png")
    for c, pg in activity_pgs.items():
        pg.to_csv(outdir / f"05_activity_gls_{c}.csv", index=False)

    candidate_periods = top_det["period_d"].head(12).to_numpy(float)
    model_summary, residuals = candidate_model_grid(det, candidate_periods, outdir / "06_sinusoid_model_fit_summary.csv")
    residuals.to_csv(outdir / "06_best_model_residuals.csv", index=False)

    lims = None
    if not args.skip_injections:
        if args.period_grid:
            periods = np.array([float(x) for x in args.period_grid.split(",")])
        else:
            periods = np.unique(np.round(np.logspace(np.log10(2.0), np.log10(min(1000.0, max_period)), 32), 6))

        if args.mass_grid:
            masses = np.array([float(x) for x in args.mass_grid.split(",")])
        else:
            masses = np.unique(np.round(np.logspace(np.log10(0.15), np.log10(20.0), 28), 6))

        part_path = run_injection_grid(
            det, outdir, periods, masses, args.n_phase,
            args.job_index, args.n_jobs, args.seed, args.mstar
        )

        if args.n_jobs == 1:
            mat = merge_injection_parts(outdir)
            lims = detection_limits_from_matrix(mat, outdir)
            plot_detection_maps(mat, lims, outdir)

    write_summary(outdir, det, det_summ, top_raw, top_det, audit, model_summary, lims)
    print("DONE")


if __name__ == "__main__":
    main()
