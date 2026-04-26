"""
Test 2 V2: MAP con priors no-circulares + componente anillo a 400d

Modelo:
  RV(t) = Kepleriano(P_p, K, e in [0,0.5], omega libre, t0)
        + Anillo_sinusoidal(P_ring~400d, A_400, phi_400)
        + GP_QP(sigma, rho, P_rot~38.65, tau)
        + offset_per_inst + jitter_per_inst

Output:
  data_interim/qc/test2_map_v2_results.json
"""
import sys, json
import numpy as np
import pandas as pd
from scipy.optimize import minimize, differential_evolution
import celerite2
from celerite2 import terms

INPUT  = "data_interim/qc/master_v3_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"

# ============================================================
# Cargar datos
# ============================================================
def load_data():
    df = pd.read_csv(INPUT)
    df = df.sort_values('bjd').reset_index(drop=True)
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    counts = df.groupby('sub_inst').size()
    keep = counts[counts >= 10].index.tolist()
    df = df[df['sub_inst'].isin(keep)].copy()
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    df['inst_idx'] = df['sub_inst'].map(
        {s: i for i, s in enumerate(sorted(df['sub_inst'].unique()))})
    return df, sorted(df['sub_inst'].unique())

# ============================================================
# Modelos
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

def build_gp_kernel(sigma, period_rot, tau):
    Q = 1.0 / np.sqrt(2)
    w0 = 2*np.pi/period_rot
    S0 = sigma**2 / (w0 * Q)
    return terms.SHOTerm(S0=S0, w0=w0, Q=Q)

# ============================================================
# Priors y posterior
# ============================================================
# Param vector layout:
# [0]  K
# [1]  P_planet
# [2]  e
# [3]  omega
# [4]  t0
# [5]  A_ring
# [6]  P_ring
# [7]  phi_ring
# [8]  log_sigma_gp
# [9]  log_rho_gp
# [10] log_period_rot
# [11] log_tau_gp
# [12..12+n_inst-1]      offsets
# [12+n_inst..]           log_jitters

def n_params(n_inst):
    return 12 + 2*n_inst

def log_posterior(params, t, rv, rv_err, inst_idx, n_inst,
                  use_planet=True, use_ring=True):
    try:
        K       = params[0]
        P_p     = params[1]
        e       = params[2]
        omega   = params[3]
        t0      = params[4]
        A_ring  = params[5]
        P_ring  = params[6]
        phi_ring= params[7]
        sigma_gp = np.exp(params[8])
        rho_gp   = np.exp(params[9])
        prot     = np.exp(params[10])
        tau_gp   = np.exp(params[11])
        offsets  = params[12:12+n_inst]
        jitters  = np.exp(params[12+n_inst:12+2*n_inst])
        
        # ===== PRIORS =====
        # K: uniforme [0, 20]
        if K < 0 or K > 20: return -np.inf
        
        # P_planet: gaussiano 232 +/- 5
        if P_p < 100 or P_p > 500: return -np.inf
        lp = -0.5 * ((P_p - 232.0)/5.0)**2
        
        # e: uniforme [0, 0.5]
        if e < 0 or e >= 0.5: return -np.inf
        
        # omega: uniforme [-pi, pi]
        if abs(omega) > np.pi: return -np.inf
        
        # t0: dentro del baseline
        if t0 < 2453000 or t0 > 2461000: return -np.inf
        
        # A_ring: uniforme [0, 10]
        if A_ring < 0 or A_ring > 10: return -np.inf
        
        # P_ring: gaussiano 400 +/- 15
        if P_ring < 350 or P_ring > 460: return -np.inf
        lp += -0.5 * ((P_ring - 400.0)/15.0)**2
        
        # phi_ring: uniforme [-pi, pi]
        if abs(phi_ring) > np.pi: return -np.inf
        
        # GP priors
        if params[8] < -3 or params[8] > 4: return -np.inf
        if params[9] < 0 or params[9] > 7: return -np.inf
        # P_rot prior estrecho
        lp += -0.5 * ((params[10] - np.log(38.65))/0.05)**2
        if params[10] < 3.0 or params[10] > 4.5: return -np.inf
        if params[11] < 0 or params[11] > 7: return -np.inf
        
        # offsets prior
        for off in offsets:
            if abs(off) > 30: return -np.inf
            lp += -0.5 * (off/15)**2
        
        # jitters
        for lj in params[12+n_inst:12+2*n_inst]:
            if lj < -3 or lj > 3: return -np.inf
        
        # ===== MODELO =====
        rv_model = np.zeros_like(t)
        if use_planet:
            rv_model += keplerian_rv(t, K, P_p, e, omega, t0)
        if use_ring:
            rv_model += ring_sinusoid(t, A_ring, P_ring, phi_ring)
        rv_model += offsets[inst_idx]
        
        resid = rv - rv_model
        yerr_total = np.sqrt(rv_err**2 + jitters[inst_idx]**2)
        
        kernel = build_gp_kernel(sigma_gp, prot, tau_gp)
        gp = celerite2.GaussianProcess(kernel, mean=0.0)
        gp.compute(t, yerr=yerr_total, quiet=True)
        ll = gp.log_likelihood(resid)
        
        if not np.isfinite(ll): return -np.inf
        return lp + ll
    except Exception:
        return -np.inf

# ============================================================
# MAP optimization con multiple restarts
# ============================================================
def run_map(t, rv, rv_err, inst_idx, n_inst, sub_inst_names,
            use_planet=True, use_ring=True, label=""):
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
    p0_base[0]  = 1.0          # K
    p0_base[1]  = 232.0        # P_planet
    p0_base[2]  = 0.2          # e (centrado en tu prior fisico)
    p0_base[3]  = 0.0          # omega
    p0_base[4]  = 2456000.0    # t0
    p0_base[5]  = 1.0          # A_ring
    p0_base[6]  = 400.0        # P_ring
    p0_base[7]  = 0.0          # phi_ring
    p0_base[8]  = np.log(2.0)  # sigma_gp
    p0_base[9]  = np.log(50.0) # rho
    p0_base[10] = np.log(38.65)# Prot
    p0_base[11] = np.log(50.0) # tau
    # offsets y jitters por defecto
    
    best_p, best_nll = p0_base.copy(), nll(p0_base)
    
    # Multiple restarts con perturbaciones
    n_restarts = 8
    for r in range(n_restarts):
        p_init = p0_base.copy()
        if r > 0:
            np.random.seed(42 + r)
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
                  f"A_ring={res.x[5]:.3f}  P_ring={res.x[6]:.1f}")
            if res.fun < best_nll:
                best_nll, best_p = res.fun, res.x
        except Exception as exc:
            print(f"  Restart {r+1}: error - {exc}")
    
    return best_p, best_nll

# ============================================================
# Main
# ============================================================
df, sub_inst_names = load_data()
n_inst = len(sub_inst_names)
t = df['bjd'].values
rv = df['rv_resid_c'].values
rv_err = np.maximum(df['rv_err'].values, 0.1)
inst_idx = df['inst_idx'].values

print(f"Data: N={len(t)}, n_inst={n_inst}")
print(f"Sub-instruments: {sub_inst_names}")

# ============================================================
# Cuatro modelos para evidencia bayesiana
# ============================================================
print(f"\n{'#'*72}")
print(f"# Comparacion de evidencia bayesiana")
print(f"#   M0: solo GP + offsets (modelo nulo)")
print(f"#   M1: GP + ring (anillo solo)")
print(f"#   M2: GP + planet (Kepleriano solo)")
print(f"#   M3: GP + planet + ring (modelo completo)")
print(f"{'#'*72}")

# M0: nulo (sin Kepleriano ni anillo)
p_M0, nll_M0 = run_map(t, rv, rv_err, inst_idx, n_inst, sub_inst_names,
                       use_planet=False, use_ring=False, label="M0 (null)")
# M1: solo anillo
p_M1, nll_M1 = run_map(t, rv, rv_err, inst_idx, n_inst, sub_inst_names,
                       use_planet=False, use_ring=True, label="M1 (ring only)")
# M2: solo Kepleriano
p_M2, nll_M2 = run_map(t, rv, rv_err, inst_idx, n_inst, sub_inst_names,
                       use_planet=True, use_ring=False, label="M2 (Kep only)")
# M3: ambos
p_M3, nll_M3 = run_map(t, rv, rv_err, inst_idx, n_inst, sub_inst_names,
                       use_planet=True, use_ring=True, label="M3 (Kep + Ring)")

# ============================================================
# Comparacion de evidencia
# ============================================================
print(f"\n{'='*72}")
print(f"COMPARACION DE EVIDENCIA BAYESIANA APROXIMADA")
print(f"{'='*72}")
print(f"{'Modelo':<25} {'-log_post':>12} {'Δnll vs M0':>12} {'log10(BF)':>10}")
print(f"{'-'*72}")
print(f"{'M0 (null)':<25} {nll_M0:>12.3f} {'reference':>12} {'reference':>10}")
print(f"{'M1 (ring only)':<25} {nll_M1:>12.3f} {nll_M0-nll_M1:>+12.3f} {(nll_M0-nll_M1)/np.log(10):>+10.3f}")
print(f"{'M2 (Kep only)':<25} {nll_M2:>12.3f} {nll_M0-nll_M2:>+12.3f} {(nll_M0-nll_M2)/np.log(10):>+10.3f}")
print(f"{'M3 (Kep + Ring)':<25} {nll_M3:>12.3f} {nll_M0-nll_M3:>+12.3f} {(nll_M0-nll_M3)/np.log(10):>+10.3f}")

# Comparaciones criticas
print(f"\nComparaciones criticas:")
print(f"  M3 vs M2 (efecto del anillo):     Δnll = {nll_M2-nll_M3:+.3f}  log10(BF) = {(nll_M2-nll_M3)/np.log(10):+.3f}")
print(f"  M3 vs M1 (efecto del Kepleriano): Δnll = {nll_M1-nll_M3:+.3f}  log10(BF) = {(nll_M1-nll_M3)/np.log(10):+.3f}")
print(f"  M2 vs M0 (Kepleriano solo):       Δnll = {nll_M0-nll_M2:+.3f}  log10(BF) = {(nll_M0-nll_M2)/np.log(10):+.3f}")
print(f"  M1 vs M0 (Anillo solo):           Δnll = {nll_M0-nll_M1:+.3f}  log10(BF) = {(nll_M0-nll_M1)/np.log(10):+.3f}")

# Reporte modelo M3 (completo)
K_M3 = p_M3[0]; P_p_M3 = p_M3[1]; e_M3 = p_M3[2]; omega_M3 = p_M3[3]
A_ring_M3 = p_M3[5]; P_ring_M3 = p_M3[6]; sigma_gp_M3 = np.exp(p_M3[8])

M_STAR = 0.972
M_p_earth = K_M3 / 0.6395 * M_STAR**(2/3) * (P_p_M3/365.25)**(1/3) * np.sqrt(1-e_M3**2)

print(f"\n{'='*72}")
print(f"PARAMETROS MAP del MODELO COMPLETO (M3)")
print(f"{'='*72}")
print(f"  Kepleriano:")
print(f"    K       = {K_M3:.3f} m/s")
print(f"    P       = {P_p_M3:.3f} d")
print(f"    e       = {e_M3:.3f}")
print(f"    omega   = {omega_M3:.3f} rad ({np.degrees(omega_M3):.1f} deg)")
print(f"    M·sini  = {M_p_earth:.3f} M_Earth")
print(f"  Anillo:")
print(f"    A       = {A_ring_M3:.3f} m/s")
print(f"    P_ring  = {P_ring_M3:.3f} d")
print(f"  GP actividad:")
print(f"    sigma   = {sigma_gp_M3:.3f} m/s")

# Save
out = {
    'M0_nll': float(nll_M0), 'M1_nll': float(nll_M1),
    'M2_nll': float(nll_M2), 'M3_nll': float(nll_M3),
    'M3_params': p_M3.tolist(),
    'M3_K': float(K_M3), 'M3_P': float(P_p_M3), 'M3_e': float(e_M3),
    'M3_omega': float(omega_M3), 'M3_M_earth': float(M_p_earth),
    'M3_A_ring': float(A_ring_M3), 'M3_P_ring': float(P_ring_M3),
    'log10BF_M3_vs_M2': float((nll_M2-nll_M3)/np.log(10)),
    'log10BF_M3_vs_M1': float((nll_M1-nll_M3)/np.log(10)),
    'log10BF_M2_vs_M0': float((nll_M0-nll_M2)/np.log(10)),
    'log10BF_M1_vs_M0': float((nll_M0-nll_M1)/np.log(10)),
    'sub_inst_names': sub_inst_names,
}
with open(f"{OUTDIR}/test2_map_v2_results.json", 'w') as f:
    json.dump(out, f, indent=2)
print(f"\n[OK] {OUTDIR}/test2_map_v2_results.json")
