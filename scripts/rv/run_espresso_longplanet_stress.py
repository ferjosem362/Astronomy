from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize

ROOT = Path(__file__).resolve().parents[2]
IN = ROOT / "data_work/rv_espresso_sflA/espresso_sflA_model_input.csv"
OUT = ROOT / "results/espresso_longplanet_stress"
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(IN)

# Flexible column detection
bjd_col = "bjd"
rv_col = "rv_ms" if "rv_ms" in df.columns else "rv_rel_ms"
err_col = "err_ms" if "err_ms" in df.columns else ("erv_ms" if "erv_ms" in df.columns else None)

df = df[[bjd_col, rv_col, err_col]].dropna().copy()
df = df.rename(columns={bjd_col:"bjd", rv_col:"rv", err_col:"err"})

t = df["bjd"].to_numpy(float)
y = df["rv"].to_numpy(float)
e0 = df["err"].to_numpy(float)

# Robust centering
t0 = np.median(t)
x = (t - t0) / 365.25
y = y - np.median(y)

def kep_rv(t, P, K, ecc, omega, M0):
    # Simple Keplerian RV, m/s.
    n = 2*np.pi/P
    M = (M0 + n*(t - t0)) % (2*np.pi)

    E = M.copy()
    for _ in range(40):
        E -= (E - ecc*np.sin(E) - M) / (1 - ecc*np.cos(E))

    nu = 2*np.arctan2(np.sqrt(1+ecc)*np.sin(E/2), np.sqrt(1-ecc)*np.cos(E/2))
    return K * (np.cos(nu + omega) + ecc*np.cos(omega))

def trend_model(coeff, x):
    # coeff low-to-high polynomial
    return sum(c * x**i for i, c in enumerate(coeff))

def sinusoid(t, P, A, phi):
    return A * np.sin(2*np.pi*(t-t0)/P + phi)

def bic(chi2, k, n):
    return chi2 + k*np.log(n)

def fit_model(deg, Plo, Phi, Kmax, emax, include_rot=True, include_short=False, err_floor=7.0):
    err = np.sqrt(e0**2 + err_floor**2)
    n = len(y)

    # params:
    # poly coeffs deg+1
    # optional rot37: Arot, phirot
    # long kep: P,K,e,omega,M0
    # optional short: P,K,e,omega,M0

    bounds = []
    for _ in range(deg+1):
        bounds.append((-500, 500))

    if include_rot:
        bounds += [(0, 20), (0, 2*np.pi)]

    bounds += [(Plo, Phi), (0, Kmax), (0, emax), (-np.pi, np.pi), (0, 2*np.pi)]

    if include_short:
        bounds += [(52, 56), (0, 3.0), (0, emax), (-np.pi, np.pi), (0, 2*np.pi)]

    def unpack(p):
        i = 0
        poly = p[i:i+deg+1]; i += deg+1
        model = trend_model(poly, x)

        rot = None
        if include_rot:
            Arot, phirot = p[i], p[i+1]; i += 2
            model += sinusoid(t, 37.0, Arot, phirot)
            rot = (Arot, phirot)

        P,K,ec,om,M0 = p[i:i+5]; i += 5
        model += kep_rv(t, P, K, ec, om, M0)
        long = (P,K,ec,om,M0)

        short = None
        if include_short:
            P2,K2,ec2,om2,M02 = p[i:i+5]; i += 5
            model += kep_rv(t, P2, K2, ec2, om2, M02)
            short = (P2,K2,ec2,om2,M02)

        return model, long, short, rot

    def obj(p):
        m, *_ = unpack(p)
        r = (y-m)/err
        return float(np.sum(r*r))

    result_de = differential_evolution(
        obj, bounds, tol=1e-8, polish=False,
        seed=42, workers=1, updating="immediate",
        maxiter=1000, popsize=20
    )
    result = minimize(obj, result_de.x, method="Nelder-Mead",
                      options={"maxiter": 20000, "xatol":1e-8, "fatol":1e-8})

    p = result.x
    m, long, short, rot = unpack(p)
    chi2 = obj(p)
    k = len(p)
    res = y - m

    row = {
        "deg": deg,
        "Plo": Plo,
        "Phi": Phi,
        "Kmax": Kmax,
        "emax": emax,
        "include_rot37": include_rot,
        "include_short54": include_short,
        "chi2": chi2,
        "bic": bic(chi2, k, n),
        "aic": chi2 + 2*k,
        "rms_ms": float(np.std(res, ddof=1)),
        "wrms_ms": float(np.sqrt(np.average(res**2, weights=1/err**2))),
        "P_long": long[0],
        "K_long": long[1],
        "e_long": long[2],
        "omega_long_deg": np.degrees(long[3]),
        "M0_long": long[4],
        "K_hits_cap": abs(long[1]-Kmax) < 1e-3,
        "e_hits_cap": abs(long[2]-emax) < 1e-3,
    }

    if short is not None:
        row.update({
            "P_short": short[0],
            "K_short": short[1],
            "e_short": short[2],
            "omega_short_deg": np.degrees(short[3]),
            "K_short_hits_cap": abs(short[1]-3.0) < 1e-3,
        })

    return row

rows = []

windows = [
    ("228_236", 224, 236),
    ("238_250", 238, 250),
    ("245_258", 245, 258),
    ("265_275", 265, 275),
    ("220_280", 220, 280),
]

for deg in [3,4]:
    for Kmax in [1.5, 2.0, 3.0, 5.0, 10.0]:
        for emax in [0.3, 0.5, 0.65, 0.8]:
            for label, lo, hi in windows:
                print(f"fit deg={deg} {label} Kmax={Kmax} emax={emax}", flush=True)
                row = fit_model(deg, lo, hi, Kmax, emax, include_rot=True, include_short=False)
                row["window"] = label
                rows.append(row)

res = pd.DataFrame(rows).sort_values("bic")
res.to_csv(OUT / "longplanet_stress_results.csv", index=False)

print()
print("=== BEST 40 ===")
print(res.head(40).to_string(index=False))
print()
print("Saved:", OUT / "longplanet_stress_results.csv")
