#!/usr/bin/env python3
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from astropy.timeseries import LombScargle


# ============================================================
# CONFIG
# ============================================================

INPUT_CANDIDATES = [
    "data_interim/qc/espresso_eso_ccf_curated_hr.csv",
    "data_processed/espresso/espresso_eso_ccf_curated_rv.csv",
    "data_processed/espresso/espresso_ccf_header_rv.csv",
]

OUTDIR = Path("results/espresso_ccf_detrend_fit_prewhiten")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_QC = Path("data_interim/qc")
OUT_QC.mkdir(parents=True, exist_ok=True)

OUT_PERIOD = OUT_QC / "periodogram_inputs"
OUT_PERIOD.mkdir(parents=True, exist_ok=True)

CANDIDATE_PERIODS = [
    0.970108077359177,
    1.028872682954664,
    1.000000000000000,
    16.300909418729493,
    18.980457431571853,
    26.657843977983536,
    32.3314560314445,
    35.67248772338015,
    46.86003032340049,
    55.267828185286334,
    89.4351158347271,
    123.38573291495129,
    170.44140877219007,
]

# Modelos multi-periodo para probar tras mirar el GLS.
MULTI_MODELS = {
    "P35": [35.67248772338015],
    "P097": [0.970108077359177],
    "P1029": [1.028872682954664],
    "P35_P16": [35.67248772338015, 16.300909418729493],
    "P35_P19": [35.67248772338015, 18.980457431571853],
    "P35_P32": [35.67248772338015, 32.3314560314445],
    "P35_P47": [35.67248772338015, 46.86003032340049],
    "P097_P35": [0.970108077359177, 35.67248772338015],
    "P1029_P35": [1.028872682954664, 35.67248772338015],
    "P35_P16_P19": [35.67248772338015, 16.300909418729493, 18.980457431571853],
    "P35_P32_P47": [35.67248772338015, 32.3314560314445, 46.86003032340049],
}

MIN_PERIOD_D = 0.5
MAX_PERIOD_FACTOR = 0.90
SAMPLES_PER_PEAK = 12
NYQUIST_FACTOR = 8

# Detrending polinómico por segmento.
# El bloque 2023PLUS tiene deriva grande: grado 2.
# El early corto: grado 1.
# Otros segmentos: grado 0 por seguridad.
DETREND_DEGREES_BY_SEGMENT = {
    "ESPRESSO_HR_2023PLUS": 2,
    "ESPRESSO_HR_2023_EARLY": 1,
    "ESPRESSO_UHR": 0,
    "ESPRESSO_OTHER": 0,
}

DEFAULT_DETREND_DEGREE = 1

# Para no dejar que outliers dominen el trend.
SIGMA_CLIP_FOR_TREND = 5.0
MAX_CLIP_ITERS = 5


# ============================================================
# UTILS
# ============================================================

def find_existing_input():
    for p in INPUT_CANDIDATES:
        pp = Path(p)
        if pp.exists():
            return pp
    raise FileNotFoundError(
        "No encuentro ningún input. Probé:\n" +
        "\n".join(f"  {p}" for p in INPUT_CANDIDATES)
    )


def pick_col(df, names, required=True):
    for n in names:
        if n in df.columns:
            return n
    if required:
        raise KeyError(f"No encuentro columna entre {names}. Columnas: {list(df.columns)}")
    return None


def robust_numeric(s):
    return pd.to_numeric(s, errors="coerce").astype(float)


def weighted_linear_fit(X, y, yerr):
    """
    Ajuste lineal ponderado.
    Devuelve beta, yhat, residuals, chi2, cov.
    """
    y = np.asarray(y, float)
    yerr = np.asarray(yerr, float)
    X = np.asarray(X, float)

    good = (
        np.isfinite(y)
        & np.isfinite(yerr)
        & (yerr > 0)
        & np.all(np.isfinite(X), axis=1)
    )

    Xg = X[good]
    yg = y[good]
    eg = yerr[good]

    if len(yg) <= Xg.shape[1]:
        raise ValueError("No hay suficientes puntos para el ajuste.")

    w = 1.0 / np.square(eg)
    sw = np.sqrt(w)

    Xw = Xg * sw[:, None]
    yw = yg * sw

    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    yhat_g = Xg @ beta
    resid_g = yg - yhat_g
    chi2 = float(np.sum((resid_g / eg) ** 2))

    # Covarianza aproximada.
    try:
        xtwx_inv = np.linalg.inv(Xg.T @ (w[:, None] * Xg))
    except np.linalg.LinAlgError:
        xtwx_inv = np.full((Xg.shape[1], Xg.shape[1]), np.nan)

    yhat = np.full(len(y), np.nan)
    resid = np.full(len(y), np.nan)
    yhat[good] = yhat_g
    resid[good] = resid_g

    return beta, yhat, resid, chi2, xtwx_inv, good


def bic(chi2, n, k):
    return float(chi2 + k * np.log(n))


def design_matrix_offsets(t, model_inst, periods):
    """
    Diseño:
    - offset por model_inst
    - sin/cos por periodo, globales
    """
    t = np.asarray(t, float)
    inst = pd.Series(model_inst).astype(str).to_numpy()
    inst_levels = sorted(pd.unique(inst))

    cols = []
    names = []

    for lev in inst_levels:
        cols.append((inst == lev).astype(float))
        names.append(f"offset[{lev}]")

    t0 = np.nanmedian(t)
    tt = t - t0

    for p in periods:
        phase = 2.0 * np.pi * tt / float(p)
        cols.append(np.sin(phase))
        names.append(f"sin[P={p:.9g}]")
        cols.append(np.cos(phase))
        names.append(f"cos[P={p:.9g}]")

    X = np.vstack(cols).T
    return X, names, t0


def fit_period_model(df, periods, y_col="rv_detrended", err_col="rv_err"):
    t = df["time"].to_numpy(float)
    y = df[y_col].to_numpy(float)
    e = df[err_col].to_numpy(float)
    model_inst = df["model_inst"].astype(str).to_numpy()

    X, names, t0 = design_matrix_offsets(t, model_inst, periods)
    beta, yhat, resid, chi2, cov, good = weighted_linear_fit(X, y, e)

    n = int(np.sum(good))
    k = int(X.shape[1])
    this_bic = bic(chi2, n, k)
    rms = float(np.nanstd(resid[good], ddof=1))
    wrms = float(np.sqrt(np.average(resid[good] ** 2, weights=1.0 / e[good] ** 2)))

    result = {
        "periods": ",".join(f"{p:.12g}" for p in periods) if periods else "none",
        "n": n,
        "k": k,
        "chi2": chi2,
        "bic": this_bic,
        "rms_mps": rms,
        "wrms_mps": wrms,
        "t0": t0,
    }

    # Amplitudes por periodo.
    param = dict(zip(names, beta))
    for p in periods:
        sname = f"sin[P={p:.9g}]"
        cname = f"cos[P={p:.9g}]"
        a = float(param[sname])
        b = float(param[cname])
        amp = float(np.sqrt(a * a + b * b))
        phase_rad = float(np.arctan2(b, a))
        result[f"amp_P_{p:.9g}_mps"] = amp
        result[f"sin_P_{p:.9g}"] = a
        result[f"cos_P_{p:.9g}"] = b
        result[f"phase_P_{p:.9g}_rad"] = phase_rad

    return result, yhat, resid, beta, names


def sigma_clip_mask(y, sigma=5.0, max_iters=5):
    y = np.asarray(y, float)
    mask = np.isfinite(y)

    for _ in range(max_iters):
        yy = y[mask]
        if len(yy) < 10:
            break
        med = np.nanmedian(yy)
        mad = np.nanmedian(np.abs(yy - med))
        scale = 1.4826 * mad if mad > 0 else np.nanstd(yy)
        if not np.isfinite(scale) or scale <= 0:
            break
        new_mask = mask & (np.abs(y - med) <= sigma * scale)
        if new_mask.sum() == mask.sum():
            break
        mask = new_mask

    return mask


def detrend_by_segment(df):
    """
    Produce:
    - rv_segment_centered
    - trend
    - rv_detrended
    - trend_degree
    """
    out = df.copy()
    out["rv_segment_centered"] = np.nan
    out["trend"] = np.nan
    out["rv_detrended"] = np.nan
    out["trend_degree"] = np.nan

    rows = []

    for inst, g in out.groupby("model_inst", dropna=False):
        idx = g.index
        t = g["time"].to_numpy(float)
        rv = g["rv"].to_numpy(float)
        err = g["rv_err"].to_numpy(float)

        good = np.isfinite(t) & np.isfinite(rv) & np.isfinite(err) & (err > 0)
        if good.sum() < 3:
            continue

        degree = DETREND_DEGREES_BY_SEGMENT.get(str(inst), DEFAULT_DETREND_DEGREE)
        degree = min(int(degree), max(0, good.sum() - 2))

        med = np.nanmedian(rv[good])
        centered = rv - med

        tt = t - np.nanmedian(t[good])
        baseline = np.nanmax(t[good]) - np.nanmin(t[good])
        if baseline > 0:
            x = tt / baseline
        else:
            x = tt

        # Fit trend on clipped centered data.
        clip = good.copy()
        if degree > 0:
            for _ in range(MAX_CLIP_ITERS):
                X_tmp = np.vstack([x ** d for d in range(degree, -1, -1)]).T
                try:
                    beta_tmp, yhat_tmp, resid_tmp, *_ = weighted_linear_fit(
                        X_tmp[clip], centered[clip], err[clip]
                    )
                except Exception:
                    break

                new_submask = sigma_clip_mask(
                    resid_tmp[clip],
                    sigma=SIGMA_CLIP_FOR_TREND,
                    max_iters=1,
                )
                clip_indices = np.where(clip)[0]
                new_clip = clip.copy()
                new_clip[clip_indices] = new_submask

                if new_clip.sum() == clip.sum():
                    clip = new_clip
                    break
                clip = new_clip

        X = np.vstack([x ** d for d in range(degree, -1, -1)]).T

        try:
            beta, yhat, resid, chi2, cov, used = weighted_linear_fit(
                X[clip], centered[clip], err[clip]
            )
            trend = X @ beta
        except Exception:
            degree = 0
            trend = np.zeros_like(centered)
            resid = centered.copy()
            chi2 = np.nan
            used = good

        detrended = centered - trend

        out.loc[idx, "rv_segment_centered"] = centered
        out.loc[idx, "trend"] = trend
        out.loc[idx, "rv_detrended"] = detrended
        out.loc[idx, "trend_degree"] = degree

        rows.append({
            "model_inst": inst,
            "n": int(good.sum()),
            "n_used_for_trend": int(np.sum(clip)),
            "degree": int(degree),
            "time_min": float(np.nanmin(t[good])),
            "time_max": float(np.nanmax(t[good])),
            "baseline_d": float(np.nanmax(t[good]) - np.nanmin(t[good])),
            "rv_median": float(med),
            "rv_raw_std": float(np.nanstd(rv[good], ddof=1)),
            "rv_segment_centered_std": float(np.nanstd(centered[good], ddof=1)),
            "rv_detrended_std": float(np.nanstd(detrended[good], ddof=1)),
            "rv_err_median": float(np.nanmedian(err[good])),
            "trend_chi2": chi2,
        })

    summary = pd.DataFrame(rows)
    return out, summary


def compute_gls(df, y_col, err_col, label, min_period=MIN_PERIOD_D):
    t = df["time"].to_numpy(float)
    y = df[y_col].to_numpy(float)
    e = df[err_col].to_numpy(float)

    good = np.isfinite(t) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    t = t[good]
    y = y[good]
    e = e[good]

    baseline = float(np.nanmax(t) - np.nanmin(t))
    max_period = max(min_period * 1.5, baseline * MAX_PERIOD_FACTOR)

    ls = LombScargle(t, y, e, center_data=True, fit_mean=True)
    freq, power = ls.autopower(
        minimum_frequency=1.0 / max_period,
        maximum_frequency=1.0 / min_period,
        samples_per_peak=SAMPLES_PER_PEAK,
        nyquist_factor=NYQUIST_FACTOR,
    )

    period = 1.0 / freq

    try:
        fap_10 = float(ls.false_alarm_level(0.10))
        fap_01 = float(ls.false_alarm_level(0.01))
        fap_001 = float(ls.false_alarm_level(0.001))
    except Exception:
        fap_10 = np.nan
        fap_01 = np.nan
        fap_001 = np.nan

    pg = pd.DataFrame({
        "period_d": period,
        "frequency_1_d": freq,
        "power": power,
        "fap_10_level": fap_10,
        "fap_01_level": fap_01,
        "fap_001_level": fap_001,
    }).sort_values("period_d")

    pg.to_csv(OUTDIR / f"{label}_gls_periodogram.csv", index=False)

    peaks = top_independent_peaks(pg, n=40)
    header = pd.DataFrame({
        "key": [
            "input_rows", "good_rows", "time_min", "time_max", "baseline_d",
            "rv_std_mps", "rv_err_median_mps", "min_period_d", "max_period_d",
            "fap_10", "fap_01", "fap_001",
        ],
        "value": [
            len(df), len(t), np.nanmin(t), np.nanmax(t), baseline,
            np.nanstd(y, ddof=1), np.nanmedian(e), min_period, max_period,
            fap_10, fap_01, fap_001,
        ]
    })
    header.to_csv(OUTDIR / f"{label}_summary.csv", index=False)
    peaks.to_csv(OUTDIR / f"{label}_top_peaks.csv", index=False)

    plot_periodogram(pg, label, OUTDIR / f"{label}_gls_periodogram.png")
    return pg, peaks


def compute_window(df, label):
    t = df["time"].to_numpy(float)
    good = np.isfinite(t)
    t = t[good]

    baseline = float(np.nanmax(t) - np.nanmin(t))
    min_period = MIN_PERIOD_D
    max_period = max(min_period * 1.5, baseline * MAX_PERIOD_FACTOR)

    # Grid manual y estrictamente creciente en frecuencia.
    fmin = 1.0 / max_period
    fmax = 1.0 / min_period
    nfreq = int(max(5000, baseline * SAMPLES_PER_PEAK * (fmax - fmin)))
    freq = np.linspace(fmin, fmax, nfreq)

    # Window = LS de una constante, sin fit_mean.
    y = np.ones_like(t)
    ls = LombScargle(t, y, fit_mean=False, center_data=False)
    power = ls.power(freq, method="slow")

    pg = pd.DataFrame({
        "period_d": 1.0 / freq,
        "frequency_1_d": freq,
        "window_power": power,
    }).sort_values("period_d")

    pg.to_csv(OUTDIR / f"{label}_window_periodogram.csv", index=False)
    peaks = top_independent_peaks(
        pg.rename(columns={"window_power": "power"}),
        n=40,
    ).rename(columns={"power": "window_power"})
    peaks.to_csv(OUTDIR / f"{label}_window_top_peaks.csv", index=False)

    plot_window(pg, label, OUTDIR / f"{label}_window_periodogram.png")
    return pg, peaks


def top_independent_peaks(pg, n=30, exclusion_frac=0.03):
    """
    Selección simple de máximos separados en periodo.
    """
    d = pg.copy()
    d = d[np.isfinite(d["period_d"]) & np.isfinite(d["power"])]
    d = d.sort_values("power", ascending=False)

    selected = []
    periods = []

    for _, r in d.iterrows():
        p = float(r["period_d"])
        if all(abs(p - q) / max(p, q) > exclusion_frac for q in periods):
            selected.append(r)
            periods.append(p)
        if len(selected) >= n:
            break

    if not selected:
        return pd.DataFrame(columns=pg.columns)

    return pd.DataFrame(selected).reset_index(drop=True)


def plot_periodogram(pg, title, outpng):
    x = pg["period_d"].to_numpy(float)
    y = pg["power"].to_numpy(float)

    plt.figure(figsize=(12, 6))
    plt.plot(x, y, lw=1)
    if "fap_10_level" in pg.columns:
        for col, lab in [
            ("fap_10_level", "FAP 10%"),
            ("fap_01_level", "FAP 1%"),
            ("fap_001_level", "FAP 0.1%"),
        ]:
            val = pg[col].dropna()
            if len(val):
                plt.axhline(float(val.iloc[0]), ls="--", lw=0.8, label=lab)

    for p in CANDIDATE_PERIODS:
        if x.min() <= p <= x.max():
            plt.axvline(p, lw=0.6, alpha=0.45)

    plt.xscale("log")
    plt.xlabel("Period [d]")
    plt.ylabel("GLS power")
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(outpng, dpi=160)
    plt.close()


def plot_window(pg, title, outpng):
    x = pg["period_d"].to_numpy(float)
    y = pg["window_power"].to_numpy(float)

    plt.figure(figsize=(12, 6))
    plt.plot(x, y, lw=1)

    for p in CANDIDATE_PERIODS:
        if x.min() <= p <= x.max():
            plt.axvline(p, lw=0.6, alpha=0.45)

    plt.xscale("log")
    plt.xlabel("Period [d]")
    plt.ylabel("Window power")
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(outpng, dpi=160)
    plt.close()


def plot_timeseries(df, y_col, outpng, title):
    plt.figure(figsize=(13, 6))
    for inst, g in df.groupby("model_inst", dropna=False):
        plt.errorbar(
            g["time"],
            g[y_col],
            yerr=g["rv_err"],
            fmt=".",
            ms=3,
            lw=0.6,
            alpha=0.8,
            label=str(inst),
        )
    plt.xlabel("BJD")
    plt.ylabel(f"{y_col} [m/s]")
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(outpng, dpi=160)
    plt.close()


def phase_plot(df, period, y_col, outpng, title):
    t = df["time"].to_numpy(float)
    y = df[y_col].to_numpy(float)
    e = df["rv_err"].to_numpy(float)
    phase = ((t - np.nanmedian(t)) / period) % 1.0

    order = np.argsort(phase)
    plt.figure(figsize=(10, 6))
    plt.errorbar(phase[order], y[order], yerr=e[order], fmt=".", ms=3, lw=0.5)
    plt.errorbar(phase[order] + 1, y[order], yerr=e[order], fmt=".", ms=3, lw=0.5)
    plt.xlabel("Phase")
    plt.ylabel(f"{y_col} [m/s]")
    plt.title(title)
    plt.grid(True, alpha=0.25)
    plt.tight_layout()
    plt.savefig(outpng, dpi=160)
    plt.close()


def local_diagnostics(gls_pg, window_pg, periods):
    rows = []

    g = gls_pg.sort_values("period_d").copy()
    w = window_pg.sort_values("period_d").copy()

    for p in periods:
        # Exact nearest
        gi = int(np.nanargmin(np.abs(g["period_d"].to_numpy(float) - p)))
        wi = int(np.nanargmin(np.abs(w["period_d"].to_numpy(float) - p)))

        exact_power = float(g.iloc[gi]["power"])
        exact_window = float(w.iloc[wi]["window_power"])

        # Local best ±3%
        lo, hi = p * 0.97, p * 1.03
        gg = g[(g["period_d"] >= lo) & (g["period_d"] <= hi)].copy()
        ww = w[(w["period_d"] >= lo) & (w["period_d"] <= hi)].copy()

        if len(gg):
            best = gg.iloc[int(np.nanargmax(gg["power"].to_numpy(float)))]
            local_p = float(best["period_d"])
            local_pow = float(best["power"])
            wj = int(np.nanargmin(np.abs(w["period_d"].to_numpy(float) - local_p)))
            local_win = float(w.iloc[wj]["window_power"])
        else:
            local_p, local_pow, local_win = np.nan, np.nan, np.nan

        rows.append({
            "target_period_d": p,
            "exact_power": exact_power,
            "exact_window_power": exact_window,
            "local_best_period_d": local_p,
            "local_best_power": local_pow,
            "local_best_window_power": local_win,
            "local_power_over_window": local_pow / local_win if local_win and np.isfinite(local_win) else np.nan,
        })

    return pd.DataFrame(rows)


# ============================================================
# LOAD / NORMALIZE
# ============================================================

inp = find_existing_input()
print("Input:", inp)

raw = pd.read_csv(inp)
print("Raw shape:", raw.shape)
print("Raw columns:", list(raw.columns))

time_col = pick_col(raw, ["time", "bjd", "bjd_approx", "BJD", "bjd_dace"])
rv_col = pick_col(raw, ["rv", "rv_mps", "RV"])
err_col = pick_col(raw, ["rv_err", "rv_err_mps", "e_rv", "RV_ERR"])
mjd_col = pick_col(raw, ["mjd_obs", "mjd", "MJD"], required=False)
inst_col = pick_col(raw, ["model_inst", "sub_inst", "inst", "instrument"], required=False)
seg_col = pick_col(raw, ["espresso_segment", "segment"], required=False)
berv_col = pick_col(raw, ["berv", "berv_mps", "cal_berv"], required=False)
fwhm_col = pick_col(raw, ["ccf_fwhm", "ccf_fwhm_mps", "ccf_fwhm_kms"], required=False)
bis_col = pick_col(raw, ["ccf_bispan", "ccf_bis_span_mps", "ccf_bispan_mps"], required=False)
contrast_col = pick_col(raw, ["ccf_contrast", "contrast"], required=False)
mask_col = pick_col(raw, ["mask", "ccf_mask"], required=False)
mode_col = pick_col(raw, ["ins_mode", "mode"], required=False)
filename_col = pick_col(raw, ["filename", "file"], required=False)

df = pd.DataFrame()
df["time"] = robust_numeric(raw[time_col])
df["rv"] = robust_numeric(raw[rv_col])
df["rv_err"] = robust_numeric(raw[err_col])

if mjd_col:
    df["mjd_obs"] = robust_numeric(raw[mjd_col])
else:
    df["mjd_obs"] = df["time"] - 2400000.5

if inst_col:
    df["model_inst"] = raw[inst_col].astype(str)
else:
    df["model_inst"] = "ESPRESSO_ALL"

if seg_col:
    df["espresso_segment"] = raw[seg_col].astype(str)
else:
    df["espresso_segment"] = df["model_inst"]

for new_col, old_col in [
    ("berv", berv_col),
    ("ccf_fwhm", fwhm_col),
    ("ccf_bispan", bis_col),
    ("ccf_contrast", contrast_col),
]:
    if old_col:
        df[new_col] = robust_numeric(raw[old_col])
    else:
        df[new_col] = np.nan

df["mask"] = raw[mask_col].astype(str) if mask_col else ""
df["ins_mode"] = raw[mode_col].astype(str) if mode_col else ""
df["filename"] = raw[filename_col].astype(str) if filename_col else ""

# Si viene de un archivo con rv_err en mm/s por error, corregimos.
# En este proyecto vimos ratio 1000 en una versión intermedia.
med_err = np.nanmedian(df["rv_err"])
if np.isfinite(med_err) and med_err > 10:
    print(f"WARNING: rv_err median={med_err:.3f}; parece estar en mm/s. Divido por 1000.")
    df["rv_err"] = df["rv_err"] / 1000.0

good = (
    np.isfinite(df["time"])
    & np.isfinite(df["rv"])
    & np.isfinite(df["rv_err"])
    & (df["rv_err"] > 0)
)
df = df[good].copy().sort_values("time").reset_index(drop=True)

print("Clean shape:", df.shape)
print("Time range:", df["time"].min(), df["time"].max())
print("RV range:", df["rv"].min(), df["rv"].max())
print("RV err median:", df["rv_err"].median())

# ============================================================
# DETREND
# ============================================================

df_det, detrend_summary = detrend_by_segment(df)

det_out = OUT_QC / "espresso_eso_ccf_curated_hr_full_detrended.csv"
per_out = OUT_PERIOD / "espresso_eso_ccf_curated_hr_full_detrended.csv"

df_det.to_csv(det_out, index=False)

df_period = df_det[[
    "time", "rv_detrended", "rv_err", "model_inst", "mjd_obs",
    "berv", "ccf_fwhm", "ccf_bispan", "ccf_contrast",
]].rename(columns={"rv_detrended": "rv"}).copy()
df_period.to_csv(per_out, index=False)

detrend_summary.to_csv(OUTDIR / "detrend_summary.csv", index=False)

plot_timeseries(
    df_det,
    "rv_segment_centered",
    OUTDIR / "timeseries_segment_centered.png",
    "ESPRESSO CCF HR - segment centered",
)

plot_timeseries(
    df_det,
    "rv_detrended",
    OUTDIR / "timeseries_detrended.png",
    "ESPRESSO CCF HR - detrended by segment",
)

print("Saved:", det_out)
print("Saved:", per_out)
print("Saved:", OUTDIR / "detrend_summary.csv")

# ============================================================
# GLS: centered, detrended, window
# ============================================================

gls_centered, peaks_centered = compute_gls(
    df_det,
    y_col="rv_segment_centered",
    err_col="rv_err",
    label="segment_centered",
)

gls_det, peaks_det = compute_gls(
    df_det,
    y_col="rv_detrended",
    err_col="rv_err",
    label="detrended",
)

window_pg, window_peaks = compute_window(
    df_det,
    label="detrended",
)

diag = local_diagnostics(gls_det, window_pg, CANDIDATE_PERIODS)
diag.to_csv(OUTDIR / "target_period_diagnostics.csv", index=False)

print("Saved:", OUTDIR / "segment_centered_gls_periodogram.csv")
print("Saved:", OUTDIR / "detrended_gls_periodogram.csv")
print("Saved:", OUTDIR / "detrended_window_periodogram.csv")
print("Saved:", OUTDIR / "target_period_diagnostics.csv")

# ============================================================
# FITS: null + single + multi
# ============================================================

fit_rows = []

null_res, null_yhat, null_resid, *_ = fit_period_model(df_det, [], y_col="rv_detrended", err_col="rv_err")
null_res["model"] = "NULL_offsets_only"
null_res["n_periods"] = 0
fit_rows.append(null_res)

# Single period models.
single_residuals = {}
for p in CANDIDATE_PERIODS:
    name = f"P_{p:.6f}d".replace(".", "p")
    res, yhat, resid, *_ = fit_period_model(df_det, [p], y_col="rv_detrended", err_col="rv_err")
    res["model"] = name
    res["n_periods"] = 1
    res["delta_bic_vs_null"] = res["bic"] - null_res["bic"]
    fit_rows.append(res)

    tmp = df_det.copy()
    tmp["model_yhat"] = yhat
    tmp["residual"] = resid
    single_residuals[name] = tmp

    # Residual GLS after removing this period.
    pg_res, pk_res = compute_gls(
        tmp.rename(columns={"residual": "rv_tmp"}),
        y_col="rv_tmp",
        err_col="rv_err",
        label=f"residual_after_{name}",
    )

    phase_plot(
        df_det,
        p,
        "rv_detrended",
        OUTDIR / f"phase_{name}.png",
        f"Phase plot detrended RV at P={p:.6f} d",
    )

# Multi-period models.
for name, periods in MULTI_MODELS.items():
    res, yhat, resid, *_ = fit_period_model(df_det, periods, y_col="rv_detrended", err_col="rv_err")
    res["model"] = name
    res["n_periods"] = len(periods)
    res["delta_bic_vs_null"] = res["bic"] - null_res["bic"]
    fit_rows.append(res)

    tmp = df_det.copy()
    tmp["model_yhat"] = yhat
    tmp["residual"] = resid
    tmp.to_csv(OUTDIR / f"residuals_after_{name}.csv", index=False)

    compute_gls(
        tmp.rename(columns={"residual": "rv_tmp"}),
        y_col="rv_tmp",
        err_col="rv_err",
        label=f"residual_after_{name}",
    )

fit_summary = pd.DataFrame(fit_rows)
fit_summary["delta_bic_vs_best"] = fit_summary["bic"] - fit_summary["bic"].min()
fit_summary = fit_summary.sort_values("bic").reset_index(drop=True)
fit_summary.to_csv(OUTDIR / "sinusoid_model_fit_summary.csv", index=False)

# Save null residuals and best model residuals.
best_model = fit_summary.iloc[0]["model"]
if best_model == "NULL_offsets_only":
    best_periods = []
elif best_model in MULTI_MODELS:
    best_periods = MULTI_MODELS[best_model]
else:
    # Extract from single models by matching approximate label.
    best_periods = []
    for p in CANDIDATE_PERIODS:
        if best_model == f"P_{p:.6f}d".replace(".", "p"):
            best_periods = [p]
            break

best_res, best_yhat, best_resid, *_ = fit_period_model(
    df_det,
    best_periods,
    y_col="rv_detrended",
    err_col="rv_err",
)

best_df = df_det.copy()
best_df["best_model"] = best_model
best_df["best_model_yhat"] = best_yhat
best_df["best_model_residual"] = best_resid
best_df.to_csv(OUTDIR / "best_model_residuals.csv", index=False)

plot_timeseries(
    best_df.rename(columns={"best_model_residual": "resid_plot"}),
    "resid_plot",
    OUTDIR / "timeseries_best_model_residuals.png",
    f"Residuals after best model: {best_model}",
)

compute_gls(
    best_df.rename(columns={"best_model_residual": "rv_tmp"}),
    y_col="rv_tmp",
    err_col="rv_err",
    label="residual_after_best_model",
)

print("Saved:", OUTDIR / "sinusoid_model_fit_summary.csv")
print("Saved:", OUTDIR / "best_model_residuals.csv")

# ============================================================
# PRINT SUMMARY
# ============================================================

print("\n==============================")
print("DETREND SUMMARY")
print("==============================")
print(detrend_summary.to_string(index=False))

print("\n==============================")
print("TOP DETRENDED GLS PEAKS")
print("==============================")
print(peaks_det.head(25).to_string(index=False))

print("\n==============================")
print("TOP WINDOW PEAKS")
print("==============================")
print(window_peaks.head(25).to_string(index=False))

print("\n==============================")
print("TARGET DIAGNOSTICS")
print("==============================")
print(diag.to_string(index=False))

print("\n==============================")
print("MODEL FIT SUMMARY")
print("==============================")
show_cols = [
    "model", "periods", "n_periods", "n", "k", "chi2", "bic",
    "delta_bic_vs_null", "delta_bic_vs_best", "rms_mps", "wrms_mps",
]
show_cols = [c for c in show_cols if c in fit_summary.columns]
print(fit_summary[show_cols].head(40).to_string(index=False))

print("\nBest model:", best_model)
print("Best periods:", best_periods)

print("\nGenerated PNGs:")
for p in sorted(OUTDIR.glob("*.png")):
    print(" ", p)

print("\nDONE")
