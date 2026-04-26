#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize
from scipy.optimize import Bounds

ROOT = Path(__file__).resolve().parents[2]
INP = ROOT / "data_work/rv_espresso_sflA/espresso_sflA_model_input.csv"
OUT = ROOT / "results/espresso_longplanet_stress_fast"
OUT.mkdir(parents=True, exist_ok=True)

PARTIAL = OUT / "longplanet_stress_fast_partial.csv"
FINAL = OUT / "longplanet_stress_fast_results.csv"

# -----------------------------
# Configuration: FAST TRACK
# -----------------------------
DEGREES = [3, 4]

WINDOWS = {
    "228_236": (224.0, 236.0),
    "238_250": (238.0, 250.0),
    "245_258": (245.0, 258.0),
    "265_275": (265.0, 275.0),
    "220_280": (220.0, 280.0),
}

KMAX_LIST = [1.5, 2.0, 3.0, 5.0]
EMAX_LIST = [0.3, 0.5, 0.65]

# Conservative error floor. Do NOT use sub-m/s internal errors for model selection here.
ERR_FLOOR_MS = 7.0

# Differential evolution fast settings.
DE_MAXITER = 120
DE_POPSIZE = 10
DE_TOL = 1e-5
SEED = 42

# -----------------------------
# Keplerian helpers
# -----------------------------
def solve_kepler(M, e, n_iter=30):
    M = np.asarray(M, dtype=float)
    E = M.copy()
    for _ in range(n_iter):
        E -= (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
    return E

def kepler_rv(t, P, K, e, omega, M0, tref):
    M = 2*np.pi*((t - tref)/P) + M0
    M = np.mod(M, 2*np.pi)
    E = solve_kepler(M, e)
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    return K * (np.cos(nu + omega) + e*np.cos(omega))

def poly_design(t, deg):
    x = (t - np.median(t)) / 365.25
    return np.vstack([x**i for i in range(deg + 1)]).T

def weighted_poly_fit(t, y, err, deg):
    X = poly_design(t, deg)
    w = 1.0 / np.asarray(err, float)
    Xw = X * w[:, None]
    yw = y * w
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    model = X @ beta
    return beta, model

def chi2_of_model(y, model, err):
    r = (y - model) / err
    return float(np.sum(r*r))

def bic_aic(chi2, n, k):
    return chi2 + k*np.log(n), chi2 + 2*k

def fit_poly_only(t, y, err, deg):
    beta, model = weighted_poly_fit(t, y, err, deg)
    resid = y - model
    chi2 = chi2_of_model(y, model, err)
    n = len(y)
    k = deg + 1
    bic, aic = bic_aic(chi2, n, k)
    return {
        "model": "poly_only",
        "deg": deg,
        "window": "none",
        "Kmax": np.nan,
        "emax": np.nan,
        "P_long": np.nan,
        "K_long": np.nan,
        "e_long": np.nan,
        "omega_deg": np.nan,
        "M0_rad": np.nan,
        "K_hits_cap": False,
        "e_hits_cap": False,
        "P_hits_edge": False,
        "chi2": chi2,
        "bic": bic,
        "aic": aic,
        "rms_ms": float(np.std(resid, ddof=1)),
        "wrms_ms": float(np.sqrt(np.average(resid**2, weights=1/err**2))),
        "n": n,
        "n_params": k,
    }, model, resid

def fit_poly_plus_kep(t, y, err, deg, window_name, Plo, Phi, Kmax, emax):
    tref = np.median(t)
    n = len(y)

    # params: P, K, e, omega, M0, poly coefficients...
    # But for speed/stability, optimize only P,K,e,omega,M0.
    # For each trial, solve weighted polynomial linearly.
    def objective(theta):
        P, K, ecc, omega, M0 = theta
        kep = kepler_rv(t, P, K, ecc, omega, M0, tref)
        beta, trend = weighted_poly_fit(t, y - kep, err, deg)
        model = trend + kep
        return chi2_of_model(y, model, err)

    bounds = [
        (Plo, Phi),
        (0.0, Kmax),
        (0.0, emax),
        (-np.pi, np.pi),
        (0.0, 2*np.pi),
    ]

    result_de = differential_evolution(
        objective,
        bounds=bounds,
        maxiter=DE_MAXITER,
        popsize=DE_POPSIZE,
        tol=DE_TOL,
        polish=False,
        seed=SEED,
        updating="immediate",
        workers=1,
    )

    result = minimize(
        objective,
        result_de.x,
        method="Nelder-Mead",
        options={"maxiter": 3000, "xatol": 1e-8, "fatol": 1e-8},
    )

    theta = result.x if result.success else result_de.x

    # Clip into legal bounds.
    P = float(np.clip(theta[0], Plo, Phi))
    K = float(np.clip(theta[1], 0.0, Kmax))
    ecc = float(np.clip(theta[2], 0.0, emax))
    omega = float(((theta[3] + np.pi) % (2*np.pi)) - np.pi)
    M0 = float(theta[4] % (2*np.pi))

    kep = kepler_rv(t, P, K, ecc, omega, M0, tref)
    beta, trend = weighted_poly_fit(t, y - kep, err, deg)
    model = trend + kep
    resid = y - model
    chi2 = chi2_of_model(y, model, err)

    k = (deg + 1) + 5
    bic, aic = bic_aic(chi2, n, k)

    epsP = 0.01 * (Phi - Plo)
    return {
        "model": "poly_plus_kep",
        "deg": deg,
        "window": window_name,
        "Kmax": Kmax,
        "emax": emax,
        "P_long": P,
        "K_long": K,
        "e_long": ecc,
        "omega_deg": float(np.degrees(omega)),
        "M0_rad": M0,
        "K_hits_cap": bool(K > 0.995*Kmax),
        "e_hits_cap": bool(ecc > 0.995*emax),
        "P_hits_edge": bool((P - Plo) < epsP or (Phi - P) < epsP),
        "chi2": chi2,
        "bic": bic,
        "aic": aic,
        "rms_ms": float(np.std(resid, ddof=1)),
        "wrms_ms": float(np.sqrt(np.average(resid**2, weights=1/err**2))),
        "n": n,
        "n_params": k,
        "success_nm": bool(result.success),
        "message_nm": str(result.message),
    }, model, resid

def save_checkpoint(rows):
    df = pd.DataFrame(rows)
    df.to_csv(PARTIAL, index=False)

    cols = [c for c in [
        "model","deg","window","Kmax","emax",
        "P_long","K_long","e_long",
        "K_hits_cap","e_hits_cap","P_hits_edge",
        "delta_bic_vs_poly","bic","aic","chi2","rms_ms","wrms_ms"
    ] if c in df.columns]

    print()
    print("="*90)
    print("CHECKPOINT rows:", len(df), "->", PARTIAL)
    if len(df):
        print(df.sort_values("bic")[cols].head(12).to_string(index=False))
    print("="*90)
    print()

def main():
    print("INPUT:", INP)
    print("OUT:", OUT)

    df = pd.read_csv(INP)
    required = ["bjd", "rv_ms"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing column {c} in {INP}")

    df = df.dropna(subset=["bjd", "rv_ms"]).copy()
    t = df["bjd"].to_numpy(float)
    y = df["rv_ms"].to_numpy(float)

    if "err_ms" in df.columns:
        err_raw = df["err_ms"].to_numpy(float)
    elif "err_ms_raw" in df.columns:
        err_raw = df["err_ms_raw"].to_numpy(float)
    else:
        err_raw = np.ones_like(y)

    err = np.sqrt(np.asarray(err_raw, float)**2 + ERR_FLOOR_MS**2)

    print("n:", len(y))
    print("BJD range:", float(np.min(t)), float(np.max(t)), "baseline:", float(np.max(t)-np.min(t)))
    print("RV mean/std/min/max:", float(np.mean(y)), float(np.std(y, ddof=1)), float(np.min(y)), float(np.max(y)))
    print("ERR median:", float(np.median(err)), "floor:", ERR_FLOOR_MS)

    rows = []
    poly_baseline = {}

    for deg in DEGREES:
        print()
        print("#"*90)
        print("DEGREE", deg)
        print("#"*90)

        row0, model0, resid0 = fit_poly_only(t, y, err, deg)
        rows.append(row0)
        poly_baseline[deg] = row0["bic"]
        save_checkpoint(rows)

        for Kmax in KMAX_LIST:
            for emax in EMAX_LIST:
                for wname, (Plo, Phi) in WINDOWS.items():
                    print(f"fit deg={deg} {wname} Kmax={Kmax} emax={emax}", flush=True)
                    row, model, resid = fit_poly_plus_kep(
                        t, y, err, deg, wname, Plo, Phi, Kmax, emax
                    )
                    row["delta_bic_vs_poly"] = poly_baseline[deg] - row["bic"]
                    rows.append(row)

                    # Save timeseries only for currently good-looking fits to avoid clutter.
                    if row["delta_bic_vs_poly"] > 5:
                        safe = f"deg{deg}_{wname}_K{Kmax}_e{emax}".replace(".", "p")
                        ts = pd.DataFrame({
                            "bjd": t,
                            "rv_ms": y,
                            "err_ms": err,
                            "model_ms": model,
                            "resid_ms": resid,
                        })
                        ts.to_csv(OUT / f"{safe}_timeseries.csv", index=False)

                    save_checkpoint(rows)

    res = pd.DataFrame(rows)
    res.to_csv(FINAL, index=False)

    print()
    print("#"*90)
    print("FINAL BEST BY BIC")
    print("#"*90)
    cols = [c for c in [
        "model","deg","window","Kmax","emax",
        "P_long","K_long","e_long",
        "K_hits_cap","e_hits_cap","P_hits_edge",
        "delta_bic_vs_poly","bic","aic","chi2","rms_ms","wrms_ms"
    ] if c in res.columns]
    print(res.sort_values("bic")[cols].head(30).to_string(index=False))
    print()
    print("SAVED:")
    print(" ", PARTIAL)
    print(" ", FINAL)

if __name__ == "__main__":
    main()
