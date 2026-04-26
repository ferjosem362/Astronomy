"""
Test 2 V3: MAP comparativo M0/M1/M2/M3 sobre dataset homogeneo HARPS+ESPRESSO.

Criterio de exclusion (declarado a priori):
  Excluir sub-instrumentos con std_resid > 5 m/s
  -> excluye CHIRON_zhao2017 (std=10.2) y ES_zhao2017 (std=7.1)

Dataset esperado: HARPS03 + HARPS15 + ESPRESSO19_HR11 + ESPRESSO_EXT_2023plus
                  ~479 noches sobre 19+ yr

Modelo y priors identicos a V2.
"""
import sys, json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import celerite2
from celerite2 import terms

INPUT  = "data_interim/qc/master_v3_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"
OUTPUT = f"{OUTDIR}/test2_map_v3_homogeneous.json"

# CRITERIO DE EXCLUSION (a priori)
STD_THRESHOLD = 5.0   # m/s
MIN_NIGHTS = 10        # noches por sub-inst

# ============================================================
def load_homogeneous_data():
    df = pd.read_csv(INPUT)
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    
    # Computar std por sub_inst
    print("\n  Filtrado homogeneo - inventario completo:")
    print(f"  {'sub_inst':<28} {'N':>5} {'std':>7}  {'incluir':<10}")
    
    inst_stats = []
    keep_subs = []
    for sub in sorted(df['sub_inst'].unique()):
        m = df['sub_inst'] == sub
        n = m.sum()
        std = df.loc[m, 'rv_resid_c'].std()
        keep = (n >= MIN_NIGHTS) and (std <= STD_THRESHOLD)
        inst_stats.append({'sub_inst': sub, 'N': n, 'std': std, 'keep': keep})
        flag = "INCLUIR" if keep else f"excluir ({'std>'+str(STD_THRESHOLD) if std>STD_THRESHOLD else 'N<'+str(MIN_NIGHTS)})"
        print(f"  {sub:<28} {n:>5} {std:>7.2f}  {flag}")
        if keep:
            keep_subs.append(sub)
    
    df = df[df['sub_inst'].isin(keep_subs)].copy()
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    df['inst_idx'] = df['sub_inst'].map({s: i for i, s in enumerate(sorted(df['sub_inst'].unique()))})
    df = df.sort_values('bjd').reset_index(drop=True)
    
    print(f"\n  Dataset homogeneo: N={len(df)}, {len(keep_subs)} sub-instrumentos")
    print(f"  Baseline: {(df['bjd'].max()-df['bjd'].min())/365.25:.2f} yr")
    return df, sorted(df['sub_inst'].unique())

# ============================================================
def kepler_solve(M, e, tol=1e-10, max_iter=200):
    E = M.copy()
    for _ in range(max_iter):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < tol: break
    return E

def keplerian_rv(t, K, P, e, omega, t0):
    M = 2*np.pi*(t - t0)/P
    E = kepler_solve(M, e)
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    return K * (np.cos(nu + omega) + e*np.cos(omega))

def ring_sinusoid(t, A, P_ring, phi):
    return A * np.sin(2*np.pi*t/P_ring + phi)

def build_gp_kernel(sigma, period_rot):
    Q = 1.0 / np.sqrt(2)
    w0 = 2*np.pi/period_rot
    S0 = sigma**2 / (w0 * Q)
    return terms.SHOTerm(S0=S0, w0=w0, Q=Q)

# ============================================================
def n_params(n_inst):
    return 12 + 2*n_inst

def log_posterior(params, t, rv, rv_err, inst_idx, n_inst,
                  use_planet=True, use_ring=True):
    try:
        K = params[0]; P_p = params[1]; e = params[2]
        omega = params[3]; t0 = params[4]
        A_ring = params[5]; P_ring = params[6]; phi_ring = params[7]
        sigma_gp = np.exp(params[8]); rho_gp = np.exp(params[9])
        prot = np.exp(params[10]); tau_gp = np.exp(params[11])
        offsets = params[12:12+n_inst]
        jitters = np.exp(params[12+n_inst:12+2*n_inst])
        
        # Priors
        if K < 0 or K > 20: return -np.inf
        if P_p < 100 or P_p > 500: return -np.inf
        lp = -0.5 * ((P_p - 232.0)/5.0)**2
        if e < 0 or e >= 0.5: return -np.inf
        if abs(omega) > np.pi: return -np.inf
        if t0 < 2453000 or t0 > 2461000: return -np.inf
        if A_ring < 0 or A_ring > 10: return -np.inf
        if P_ring < 350 or P_ring > 460: return -np.inf
        lp += -0.5 * ((P_ring - 400.0)/15.0)**2
        if abs(phi_ring) > np.pi: return -np.inf
        if params[8] < -3 or params[8] > 4: return -np.inf
        if params[9] < 0 or params[9] > 7: return -np.inf
        lp += -0.5 * ((params[10] - np.log(38.65))/0.05)**2
        if params[10] < 3.0 or params[10] > 4.5: return -np.inf
        if params[11] < 0 or params[11] > 7: return -np.inf
        for off in offsets:
            if abs(off) > 30: return -np.inf
            lp += -0.5 * (off/15)**2
        for lj in params[12+n_inst:12+2*n_inst]:
            if lj < -3 or lj > 3: return -np.inf
        
        # Modelo
        rv_model = np.zeros_like(t)
        if use_planet:
            rv_model += keplerian_rv(t, K, P_p, e, omega, t0)
        if use_ring:
            rv_model += ring_sinusoid(t, A_ring, P_ring, phi_ring)
        rv_model += offsets[inst_idx]
        
        resid = rv - rv_model
        yerr_total = np.sqrt(rv_err**2 + jitters[inst_idx]**2)
        
        kernel = build_gp_kernel(sigma_gp, prot)
        gp = celerite2.GaussianProcess(kernel, mean=0.0)
        gp.compute(t, yerr=yerr_total, quiet=True)
        ll = gp.log_likelihood(resid)
        if not np.isfinite(ll): return -np.inf
        return lp + ll
    except Exception:
        return -np.inf

# ============================================================
def run_map(t, rv, rv_err, inst_idx, n_inst, use_planet=True, use_ring=True, label=""):
    n_p = n_params(n_inst)
    
    def nll(p):
        lp = log_posterior(p, t, rv, rv_err, inst_idx, n_inst,
                          use_planet=use_planet, use_ring=use_ring)
        if not np.isfinite(lp): return 1e30
        return -lp
    
    print(f"\n{'='*72}")
    print(f"MAP run: {label}  (planet={use_planet}, ring={use_ring})")
    print(f"{'='*72}")
    
    p0_base = np.zeros(n_p)
    p0_base[0]  = 1.0
    p0_base[1]  = 232.0
    p0_base[2]  = 0.2
    p0_base[3]  = 0.0
    p0_base[4]  = 2456000.0
    p0_base[5]  = 1.0
    p0_base[6]  = 400.0
    p0_base[7]  = 0.0
    p0_base[8]  = np.log(2.0)
    p0_base[9]  = np.log(50.0)
    p0_base[10] = np.log(38.65)
    p0_base[11] = np.log(50.0)
    
    best_p, best_nll = p0_base.copy(), nll(p0_base)
    
    for r in range(8):
        p_init = p0_base.copy()
        if r > 0:
            np.random.seed(42+r)
            p_init[0] = np.random.uniform(0.3, 4.0)
            p_init[1] = np.random.normal(232, 4)
            p_init[2] = np.random.uniform(0.05, 0.45)
            p_init[3] = np.random.uniform(-np.pi, np.pi)
            p_init[4] = 2456000 + np.random.uniform(-200, 200)
            p_init[5] = np.random.uniform(0.3, 3.0)
            p_init[6] = np.random.normal(400, 10)
            p_init[7] = np.random.uniform(-np.pi, np.pi)
            p_init[8] = np.log(np.random.uniform(0.5, 5))
            p_init[9] = np.log(np.random.uniform(20, 200))
            p_init[11] = np.log(np.random.uniform(20, 200))
        
        try:
            res = minimize(nll, p_init, method='Nelder-Mead',
                          options={'maxiter': 4000, 'xatol': 1e-4, 'fatol': 1e-4})
            print(f"  Restart {r+1}: nll={res.fun:.3f}  K={res.x[0]:.3f}  "
                  f"P_p={res.x[1]:.2f}  e={res.x[2]:.3f}  "
                  f"A_r={res.x[5]:.3f}  P_r={res.x[6]:.1f}  "
                  f"sigma_gp={np.exp(res.x[8]):.2f}")
            if res.fun < best_nll:
                best_nll, best_p = res.fun, res.x
        except Exception as exc:
            print(f"  Restart {r+1}: error - {exc}")
    
    return best_p, best_nll

# ============================================================
df, sub_inst_names = load_homogeneous_data()
n_inst = len(sub_inst_names)
t = df['bjd'].values
rv = df['rv_resid_c'].values
rv_err = np.maximum(df['rv_err'].values, 0.1)
inst_idx = df['inst_idx'].values

print(f"\n{'#'*72}")
print(f"# Test 2 V3 - Comparacion bayesiana 4 modelos sobre dataset homogeneo")
print(f"{'#'*72}")

p_M0, nll_M0 = run_map(t, rv, rv_err, inst_idx, n_inst, False, False, "M0 (null)")
p_M1, nll_M1 = run_map(t, rv, rv_err, inst_idx, n_inst, False, True,  "M1 (ring only)")
p_M2, nll_M2 = run_map(t, rv, rv_err, inst_idx, n_inst, True,  False, "M2 (Kep only)")
p_M3, nll_M3 = run_map(t, rv, rv_err, inst_idx, n_inst, True,  True,  "M3 (Kep + Ring)")

print(f"\n{'='*72}")
print(f"COMPARACION DE EVIDENCIA BAYESIANA - DATASET HOMOGENEO HARPS+ESPRESSO")
print(f"{'='*72}")
print(f"{'Modelo':<25} {'-log_post':>12} {'Δnll vs M0':>12} {'log10(BF)':>10}")
print(f"{'-'*72}")
print(f"{'M0 (null)':<25} {nll_M0:>12.3f} {'reference':>12} {'reference':>10}")
print(f"{'M1 (ring only)':<25} {nll_M1:>12.3f} {nll_M0-nll_M1:>+12.3f} {(nll_M0-nll_M1)/np.log(10):>+10.3f}")
print(f"{'M2 (Kep only)':<25} {nll_M2:>12.3f} {nll_M0-nll_M2:>+12.3f} {(nll_M0-nll_M2)/np.log(10):>+10.3f}")
print(f"{'M3 (Kep + Ring)':<25} {nll_M3:>12.3f} {nll_M0-nll_M3:>+12.3f} {(nll_M0-nll_M3)/np.log(10):>+10.3f}")

print(f"\nComparaciones criticas:")
print(f"  M3 vs M2 (anillo sobre Kep):    Δ={nll_M2-nll_M3:+.3f}  log10(BF)={(nll_M2-nll_M3)/np.log(10):+.3f}")
print(f"  M3 vs M1 (Kep sobre anillo):    Δ={nll_M1-nll_M3:+.3f}  log10(BF)={(nll_M1-nll_M3)/np.log(10):+.3f}")
print(f"  M2 vs M0 (Kep solo vs nulo):    Δ={nll_M0-nll_M2:+.3f}  log10(BF)={(nll_M0-nll_M2)/np.log(10):+.3f}")
print(f"  M1 vs M0 (Anillo solo vs nulo): Δ={nll_M0-nll_M1:+.3f}  log10(BF)={(nll_M0-nll_M1)/np.log(10):+.3f}")

# Reporte modelo M3
K_M3 = p_M3[0]; P_p_M3 = p_M3[1]; e_M3 = p_M3[2]; omega_M3 = p_M3[3]
A_ring_M3 = p_M3[5]; P_ring_M3 = p_M3[6]; sigma_gp_M3 = np.exp(p_M3[8])

M_STAR = 0.972
M_p_earth = K_M3 / 0.6395 * M_STAR**(2/3) * (P_p_M3/365.25)**(1/3) * np.sqrt(1-e_M3**2)

print(f"\n{'='*72}")
print(f"PARAMETROS MAP M3 (modelo completo)")
print(f"{'='*72}")
print(f"  Kepleriano: K={K_M3:.3f} m/s  P={P_p_M3:.3f} d  e={e_M3:.3f}  omega={np.degrees(omega_M3):+.1f}deg")
print(f"  M·sini   = {M_p_earth:.3f} M_Earth")
print(f"  Anillo:    A={A_ring_M3:.3f} m/s  P={P_ring_M3:.3f} d")
print(f"  GP:        sigma={sigma_gp_M3:.3f} m/s")

# Reporte tambien M2 (Kep solo)
K_M2 = p_M2[0]; P_p_M2 = p_M2[1]; e_M2 = p_M2[2]
sigma_gp_M2 = np.exp(p_M2[8])
M_p_M2 = K_M2 / 0.6395 * M_STAR**(2/3) * (P_p_M2/365.25)**(1/3) * np.sqrt(1-e_M2**2)
print(f"\nM2 (Kep solo): K={K_M2:.3f} m/s  P={P_p_M2:.3f} d  e={e_M2:.3f}  M={M_p_M2:.3f} M_Earth  sigma_gp={sigma_gp_M2:.3f}")

out = {
    'criterio_exclusion': f'std<{STD_THRESHOLD} m/s, N>={MIN_NIGHTS}',
    'sub_inst_incluidos': sub_inst_names,
    'N_total': len(t),
    'baseline_yr': float((t.max()-t.min())/365.25),
    'M0_nll': float(nll_M0), 'M1_nll': float(nll_M1),
    'M2_nll': float(nll_M2), 'M3_nll': float(nll_M3),
    'M3_K': float(K_M3), 'M3_P': float(P_p_M3), 'M3_e': float(e_M3),
    'M3_omega': float(omega_M3), 'M3_M_earth': float(M_p_earth),
    'M3_A_ring': float(A_ring_M3), 'M3_P_ring': float(P_ring_M3),
    'M3_sigma_gp': float(sigma_gp_M3),
    'M2_K': float(K_M2), 'M2_P': float(P_p_M2), 'M2_e': float(e_M2),
    'M2_M_earth': float(M_p_M2), 'M2_sigma_gp': float(sigma_gp_M2),
    'log10BF_M3_vs_M2': float((nll_M2-nll_M3)/np.log(10)),
    'log10BF_M2_vs_M0': float((nll_M0-nll_M2)/np.log(10)),
    'log10BF_M1_vs_M0': float((nll_M0-nll_M1)/np.log(10)),
}
with open(OUTPUT, 'w') as f:
    json.dump(out, f, indent=2)
print(f"\n[OK] {OUTPUT}")
