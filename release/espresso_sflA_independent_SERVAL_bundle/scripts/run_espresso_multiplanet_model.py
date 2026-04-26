#!/usr/bin/env python3
from pathlib import Path
import json
import warnings
import numpy as np
import pandas as pd

from scipy.optimize import differential_evolution, least_squares

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(".").resolve()
INFILE = ROOT / "data_work/rv_espresso_sflA/espresso_sflA_model_input.csv"
OUTDIR = ROOT / "results/espresso_multiplanet"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Configuración prudente
# ------------------------------------------------------------
ERR_FLOOR_MS = 5.0       # muy importante: no creer errores internos de 0.1 m/s
E_MAX = 0.65             # permitimos excentricidad moderada/alta, pero no absurda
RANDOM_SEED = 42

DE_MAXITER = 350
DE_POPSIZE = 10
DE_TOL = 1e-7

# Ventanas candidatas
LONG_WINDOWS = [
    ("long_228_236", 224.0, 236.0, 5.0),
    ("long_238_250", 238.0, 250.0, 5.0),
    ("long_245_258", 245.0, 258.0, 5.0),
    ("long_265_275", 265.0, 275.0, 5.0),
    ("long_220_280", 220.0, 280.0, 5.0),
]

SHORT_WINDOWS = [
    ("short_52_56", 52.0, 56.0, 3.0),
]

DEGREES = [3, 4]


# ------------------------------------------------------------
# Utilidades
# ------------------------------------------------------------
def pick_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"No encuentro ninguna columna entre {candidates}. Columnas disponibles: {list(df.columns)}")


def load_data():
    df = pd.read_csv(INFILE)

    bjd_col = pick_col(df, ["bjd", "BJD", "time", "jd"])
    rv_col = pick_col(df, ["rv_ms", "rv_rel_ms", "RV", "rv"])
    err_col = pick_col(df, ["err_ms", "erv_ms", "e_rv_ms", "rv_err_ms", "err"])

    d = df[[bjd_col, rv_col, err_col]].copy()
    d.columns = ["bjd", "rv_ms", "err_ms"]
    d = d.replace([np.inf, -np.inf], np.nan).dropna()

    d["err_ms"] = np.sqrt(np.asarray(d["err_ms"], float)**2 + ERR_FLOOR_MS**2)
    d = d.sort_values("bjd").reset_index(drop=True)

    return d


def solve_kepler(M, e, n_iter=30):
    M = np.asarray(M)
    E = M.copy()
    for _ in range(n_iter):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        E -= f / fp
    return E


def kepler_rv(t, P, K, e, omega, M0, tref):
    """
    RV = K [cos(nu + omega) + e cos omega]
    """
    M = 2.0 * np.pi * ((t - tref) / P) + M0
    M = np.mod(M, 2.0 * np.pi)

    E = solve_kepler(M, e)
    nu = 2.0 * np.arctan2(
        np.sqrt(1.0 + e) * np.sin(E / 2.0),
        np.sqrt(1.0 - e) * np.cos(E / 2.0),
    )

    return K * (np.cos(nu + omega) + e * np.cos(omega))


def design_matrix(t, deg, include_rot37=True):
    """
    Parte lineal:
    - polinomio temporal
    - actividad/rotación 37 d como sin/cos fijo
    """
    x = (t - np.mean(t)) / (0.5 * (np.max(t) - np.min(t)))
    cols = []

    for k in range(deg + 1):
        cols.append(x**k)

    if include_rot37:
        Prot = 37.0
        phase = 2.0 * np.pi * (t - t.min()) / Prot
        cols.append(np.sin(phase))
        cols.append(np.cos(phase))

    return np.vstack(cols).T


def unpack_components(theta, components):
    out = []
    i = 0
    for comp in components:
        label, p_lo, p_hi, kmax = comp
        P = theta[i + 0]
        K = theta[i + 1]
        e = theta[i + 2]
        omega = theta[i + 3]
        M0 = theta[i + 4]
        out.append((label, P, K, e, omega, M0))
        i += 5
    return out


def make_bounds(components):
    bounds = []
    for label, p_lo, p_hi, kmax in components:
        bounds.extend([
            (p_lo, p_hi),             # P
            (0.0, kmax),              # K
            (0.0, E_MAX),             # e
            (-np.pi, np.pi),          # omega
            (0.0, 2.0 * np.pi),       # M0
        ])
    return bounds


def eval_model(theta, t, y, e, deg, components, include_rot37=True):
    tref = np.median(t)

    kep_total = np.zeros_like(y)
    for label, P, K, ecc, omega, M0 in unpack_components(theta, components):
        kep_total += kepler_rv(t, P, K, ecc, omega, M0, tref)

    A = design_matrix(t, deg, include_rot37=include_rot37)

    # Ajuste lineal ponderado de polinomio + rotación, dado el Kepleriano.
    yw = (y - kep_total) / e
    Aw = A / e[:, None]

    try:
        beta, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    except np.linalg.LinAlgError:
        beta = np.zeros(A.shape[1])

    linear = A @ beta
    model = linear + kep_total
    resid = y - model
    chi2 = float(np.sum((resid / e)**2))

    return model, resid, chi2, beta


def objective(theta, t, y, e, deg, components, include_rot37=True):
    _, _, chi2, _ = eval_model(theta, t, y, e, deg, components, include_rot37)
    if not np.isfinite(chi2):
        return 1e100
    return chi2


def fit_model(name, deg, components, t, y, e, include_rot37=True):
    bounds = make_bounds(components)

    if len(bounds) == 0:
        theta = np.array([])
        model, resid, chi2, beta = eval_model(theta, t, y, e, deg, components, include_rot37)
    else:
        result_de = differential_evolution(
            objective,
            bounds=bounds,
            args=(t, y, e, deg, components, include_rot37),
            seed=RANDOM_SEED,
            maxiter=DE_MAXITER,
            popsize=DE_POPSIZE,
            polish=False,
            tol=DE_TOL,
            workers=1,
            updating="immediate",
        )

        lb = np.array([b[0] for b in bounds], dtype=float)
        ub = np.array([b[1] for b in bounds], dtype=float)

        def residual_func(theta):
            model, resid, chi2, beta = eval_model(theta, t, y, e, deg, components, include_rot37)
            return resid / e

        result_ls = least_squares(
            residual_func,
            result_de.x,
            bounds=(lb, ub),
            max_nfev=20000,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
        )

        theta = result_ls.x
        model, resid, chi2, beta = eval_model(theta, t, y, e, deg, components, include_rot37)

    n = len(y)
    n_linear = len(beta)
    n_nonlin = len(theta)
    k = n_linear + n_nonlin
    bic = chi2 + k * np.log(n)
    aic = chi2 + 2 * k
    rms = float(np.std(resid, ddof=1))
    wrms = float(np.sqrt(np.average(resid**2, weights=1.0 / e**2)))

    rows = []
    comps = unpack_components(theta, components)
    for label, P, K, ecc, omega, M0 in comps:
        rows.append({
            "model": name,
            "deg": deg,
            "component": label,
            "P_d": P,
            "K_ms": K,
            "ecc": ecc,
            "omega_deg": np.degrees(omega),
            "M0_rad": M0,
            "K_hits_cap": bool(abs(K - [c[3] for c in components if c[0] == label][0]) < 1e-4),
        })

    summary = {
        "model": name,
        "deg": deg,
        "n": n,
        "n_params": k,
        "chi2": chi2,
        "bic": bic,
        "aic": aic,
        "rms_ms": rms,
        "wrms_ms": wrms,
        "include_rot37": include_rot37,
        "n_keplerians": len(components),
    }

    return summary, rows, model, resid, beta, theta


def save_fit_timeseries(name, deg, t, y, e, model, resid):
    safe = name.replace("/", "_")
    out = pd.DataFrame({
        "bjd": t,
        "rv_ms": y,
        "err_ms": e,
        "model_ms": model,
        "resid_ms": resid,
    })
    out.to_csv(OUTDIR / f"{safe}_deg{deg}_timeseries.csv", index=False)


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    d = load_data()
    t = d["bjd"].to_numpy(float)
    y = d["rv_ms"].to_numpy(float)
    e = d["err_ms"].to_numpy(float)

    print("INPUT:", INFILE)
    print("OUTDIR:", OUTDIR)
    print("n nights:", len(d))
    print("BJD:", t.min(), t.max(), "baseline:", t.max() - t.min())
    print("RV rms:", np.std(y, ddof=1), "m/s")
    print("ERR median with floor:", np.median(e), "m/s")
    print()

    summaries = []
    component_rows = []
    saved = {}

    for deg in DEGREES:
        print("=" * 80)
        print(f"DEGREE {deg}")

        # M0: tendencia sola
        print("Fit M0 trend only")
        s, rows, model, resid, beta, theta = fit_model(
            name="M0_trend",
            deg=deg,
            components=[],
            t=t, y=y, e=e,
            include_rot37=False,
        )
        summaries.append(s)
        save_fit_timeseries("M0_trend", deg, t, y, e, model, resid)

        # M0r: tendencia + rotación 37
        print("Fit M0r trend + rot37")
        s, rows, model, resid, beta, theta = fit_model(
            name="M0r_trend_rot37",
            deg=deg,
            components=[],
            t=t, y=y, e=e,
            include_rot37=True,
        )
        summaries.append(s)
        save_fit_timeseries("M0r_trend_rot37", deg, t, y, e, model, resid)

        # Un Kepleriano largo
        for lw in LONG_WINDOWS:
            label, lo, hi, kmax = lw
            name = f"M1_rot37_{label}"
            print(f"Fit {name}: P=[{lo},{hi}], Kmax={kmax}")
            s, rows, model, resid, beta, theta = fit_model(
                name=name,
                deg=deg,
                components=[lw],
                t=t, y=y, e=e,
                include_rot37=True,
            )
            summaries.append(s)
            component_rows.extend(rows)
            save_fit_timeseries(name, deg, t, y, e, model, resid)

        # Largo + 54 d
        for lw in LONG_WINDOWS:
            for sw in SHORT_WINDOWS:
                name = f"M2_rot37_{lw[0]}_{sw[0]}"
                print(f"Fit {name}")
                s, rows, model, resid, beta, theta = fit_model(
                    name=name,
                    deg=deg,
                    components=[lw, sw],
                    t=t, y=y, e=e,
                    include_rot37=True,
                )
                summaries.append(s)
                component_rows.extend(rows)
                save_fit_timeseries(name, deg, t, y, e, model, resid)

    summ = pd.DataFrame(summaries)
    comp = pd.DataFrame(component_rows)

    # Comparación contra M0r del mismo grado
    summ["delta_bic_vs_M0r_same_deg"] = np.nan
    summ["delta_aic_vs_M0r_same_deg"] = np.nan

    for deg in DEGREES:
        base = summ[(summ["deg"] == deg) & (summ["model"] == "M0r_trend_rot37")]
        if len(base):
            b_bic = float(base.iloc[0]["bic"])
            b_aic = float(base.iloc[0]["aic"])
            idx = summ["deg"] == deg
            summ.loc[idx, "delta_bic_vs_M0r_same_deg"] = b_bic - summ.loc[idx, "bic"]
            summ.loc[idx, "delta_aic_vs_M0r_same_deg"] = b_aic - summ.loc[idx, "aic"]

    summ = summ.sort_values(["bic", "chi2"]).reset_index(drop=True)

    summ.to_csv(OUTDIR / "model_comparison.csv", index=False)
    comp.to_csv(OUTDIR / "component_params.csv", index=False)

    print()
    print("=" * 80)
    print("BEST MODELS BY BIC")
    print(summ.head(20).to_string(index=False))

    if len(comp):
        print()
        print("=" * 80)
        print("COMPONENT PARAMETERS FOR BEST MODELS")
        best_models = set(summ.head(10)["model"])
        print(comp[comp["model"].isin(best_models)].sort_values(["model", "component"]).to_string(index=False))

    print()
    print("Saved:")
    print(" ", OUTDIR / "model_comparison.csv")
    print(" ", OUTDIR / "component_params.csv")
    print(" ", OUTDIR / "*_timeseries.csv")


if __name__ == "__main__":
    main()
