"""
Test 2 V3 FINAL: comparacion bayesiana 4 modelos con:
  - Dataset homogeneo (excluir CHIRON, ES_zhao - std > 5 m/s)
  - Excentricidad libre [0, 0.5] (no Beta(1,5))
  - Periodo Kepleriano P=232 +/- 15 d (mas libre)
  - Anillo P_ring=400-500 d uniforme (banda difusa)
  - Formula RV-masa CORRECTA via utils/rv_mass.py

Modelos:
  M0: GP + offsets (nulo)
  M1: GP + ring (anillo solo)
  M2: GP + Kepleriano excentrico (planet solo)
  M3: GP + Kepleriano + ring (modelo completo)
"""
import sys, json
sys.path.insert(0, '.')

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import celerite2
from celerite2 import terms
from utils.rv_mass import msini_from_K_earth

INPUT  = "data_interim/qc/master_v3_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"
OUTPUT = f"{OUTDIR}/test2_map_v3_final.json"

STD_THRESHOLD = 5.0
MIN_NIGHTS = 10

# ============================================================
def load_homogeneous_data():
    df = pd.read_csv(INPUT)
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    
    print("\n  Filtrado homogeneo (std<=5 m/s, N>=10):")
    print(f"  {'sub_inst':<28} {'N':>5} {'std':>7}  {'incluir':<10}")
    
    keep_subs = []
    for sub in sorted(df['sub_inst'].unique()):
        m = df['sub_inst'] == sub
        n = m.sum()
        std = df.loc[m, 'rv_resid_c'].std()
        keep = (n >= MIN_NIGHTS) and (std <= STD_THRESHOLD)
        flag = "INCLUIR" if keep else f"excluir"
        print(f"  {sub:<28} {n:>5} {std:>7.2f}  {flag}")
        if keep:
            keep_subs.append(sub)
    
    df = df[df['sub_inst'].isin(keep_subs)].copy()
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    df['inst_idx'] = df['sub_inst'].map({s: i for i, s in enumerate(sorted(df['sub_inst'].unique()))})
    df = df.sort_values('bjd').reset_index(drop=True)
    
    print(f"\n  Dataset: N={len(df)}, {len(keep_subs)} sub-inst, "
          f"baseline {(df['bjd'].max()-df['bjd'].min())/365.25:.2f} yr")
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
    """
    Priors V3 FINAL:
      P_planet: 217-247 d uniforme (232 +/- 15)
      e_planet: 0-0.5 uniforme (no Beta)
      P_ring: 400-500 d uniforme (banda difusa)
    """
    try:
        K = params[0]; P_p = params[1]; e = params[2]
        omega = params[3]; t0 = params[4]
        A_ring = params[5]; P_ring = params[6]; phi_ring = params[7]
        sigma_gp = np.exp(params[8]); rho_gp = np.exp(params[9])
        prot = np.exp(params[10]); tau_gp = np.exp(params[11])
        offsets = params[12:12+n_inst]
        jitters = np.exp(params[12+n_inst:12+2*n_inst])
        
        # === PRIORS V3 FINAL ===
        if K < 0 or K > 20: return -np.inf
        # P_planet uniforme [217, 247]
        if P_p < 217 or P_p > 247: return -np.inf
        # e uniforme [0, 0.5]
        if e < 0 or e >= 0.5: return -np.inf
        if abs(omega) > np.pi: return -np.inf
        if t0 < 2453000 or t0 > 2461000: return -np.inf
        # Anillo
        if A_ring < 0 or A_ring > 10: return -np.inf
        # P_ring uniforme [400, 500] - banda difusa
        if P_ring < 400 or P_ring > 500: return -np.inf
        if abs(phi_ring) > np.pi: return -np.inf
        # GP
        if params[8] < -3 or params[8] > 4: return -np.inf
        if params[9] < 0 or params[9] > 7: return -np.inf
        lp = -0.5 * ((params[10] - np.log(38.65))/0.05)**2
        if params[10] < 3.0 or params[10] > 4.5: return -np.inf
        if params[11] < 0 or params[11] > 7: return -np.inf
        for off in offsets:
            if abs(off) > 30: return -np.inf
            lp += -0.5 * (off/15)**2
        for lj in params[12+n_inst:12+2*n_inst]:
            if lj < -3 or lj > 3: return -np.inf
        
        # === MODELO ===
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
def run_map(t, rv, rv_err, inst_idx, n_inst, use_planet, use_ring, label):
    n_p = n_params(n_inst)
    
    def nll(p):
        lp = log_posterior(p, t, rv, rv_err, inst_idx, n_inst,
                          use_planet=use_planet, use_ring=use_ring)
        if not np.isfinite(lp): return 1e30
        return -lp
    
    print(f"\n{'='*72}\nMAP: {label}\n{'='*72}")
    
    p0 = np.zeros(n_p)
    p0[0]=1.0; p0[1]=232.0; p0[2]=0.2; p0[3]=0.0; p0[4]=2456000.0
    p0[5]=1.0; p0[6]=450.0; p0[7]=0.0  # P_ring centro 450
    p0[8]=np.log(2.0); p0[9]=np.log(50.0)
    p0[10]=np.log(38.65); p0[11]=np.log(50.0)
    
    best_p, best_nll = p0.copy(), nll(p0)
    
    n_restarts = 10  # mas restarts dada amplitud de priors
    for r in range(n_restarts):
        p_init = p0.copy()
        if r > 0:
            np.random.seed(42+r)
            p_init[0] = np.random.uniform(0.3, 4.0)
            p_init[1] = np.random.uniform(220, 245)  # uniforme en prior
            p_init[2] = np.random.uniform(0.05, 0.45)
            p_init[3] = np.random.uniform(-np.pi, np.pi)
            p_init[4] = 2456000 + np.random.uniform(-200, 200)
            p_init[5] = np.random.uniform(0.3, 3.0)
            p_init[6] = np.random.uniform(410, 490)  # uniforme prior anillo
            p_init[7] = np.random.uniform(-np.pi, np.pi)
            p_init[8] = np.log(np.random.uniform(0.5, 5))
            p_init[9] = np.log(np.random.uniform(20, 200))
            p_init[11] = np.log(np.random.uniform(20, 200))
        
        try:
            res = minimize(nll, p_init, method='Nelder-Mead',
                          options={'maxiter': 5000, 'xatol': 1e-4, 'fatol': 1e-4})
            sigma_gp_r = np.exp(res.x[8])
            print(f"  R{r+1:2d}: nll={res.fun:8.3f}  K={res.x[0]:5.3f}  "
                  f"P={res.x[1]:6.2f}  e={res.x[2]:5.3f}  "
                  f"A_r={res.x[5]:5.2f}  P_r={res.x[6]:6.1f}  "
                  f"sig_gp={sigma_gp_r:5.2f}")
            if res.fun < best_nll:
                best_nll, best_p = res.fun, res.x
        except Exception as exc:
            print(f"  R{r+1}: error - {exc}")
    
    return best_p, best_nll

# ============================================================
df, sub_inst_names = load_homogeneous_data()
n_inst = len(sub_inst_names)
t = df['bjd'].values
rv = df['rv_resid_c'].values
rv_err = np.maximum(df['rv_err'].values, 0.1)
inst_idx = df['inst_idx'].values

print(f"\n{'#'*72}")
print(f"# TEST 2 V3 FINAL")
print(f"# Priors: P=232+/-15d, e libre [0,0.5], P_ring=400-500d uniforme")
print(f"# Formula RV-masa: utils/rv_mass.py (factor 7.149 corregido)")
print(f"{'#'*72}")

p_M0, nll_M0 = run_map(t, rv, rv_err, inst_idx, n_inst, False, False, "M0 (null)")
p_M1, nll_M1 = run_map(t, rv, rv_err, inst_idx, n_inst, False, True,  "M1 (ring)")
p_M2, nll_M2 = run_map(t, rv, rv_err, inst_idx, n_inst, True,  False, "M2 (Kep)")
p_M3, nll_M3 = run_map(t, rv, rv_err, inst_idx, n_inst, True,  True,  "M3 (Kep+ring)")

# ============================================================
print(f"\n{'='*72}")
print(f"COMPARACION DE EVIDENCIA BAYESIANA")
print(f"{'='*72}")
print(f"{'Modelo':<25} {'-log_post':>12} {'Δnll vs M0':>12} {'log10(BF)':>10}")
print(f"{'-'*72}")
print(f"{'M0 (null)':<25} {nll_M0:>12.3f} {'reference':>12} {'reference':>10}")
print(f"{'M1 (ring)':<25} {nll_M1:>12.3f} {nll_M0-nll_M1:>+12.3f} {(nll_M0-nll_M1)/np.log(10):>+10.3f}")
print(f"{'M2 (Kep)':<25} {nll_M2:>12.3f} {nll_M0-nll_M2:>+12.3f} {(nll_M0-nll_M2)/np.log(10):>+10.3f}")
print(f"{'M3 (Kep+ring)':<25} {nll_M3:>12.3f} {nll_M0-nll_M3:>+12.3f} {(nll_M0-nll_M3)/np.log(10):>+10.3f}")

print(f"\nComparaciones criticas:")
print(f"  M3 vs M2 (anillo sobre Kep): Δ={nll_M2-nll_M3:+.3f}  log10(BF)={(nll_M2-nll_M3)/np.log(10):+.3f}")
print(f"  M3 vs M1 (Kep sobre anillo): Δ={nll_M1-nll_M3:+.3f}  log10(BF)={(nll_M1-nll_M3)/np.log(10):+.3f}")
print(f"  M2 vs M0 (Kep solo):         Δ={nll_M0-nll_M2:+.3f}  log10(BF)={(nll_M0-nll_M2)/np.log(10):+.3f}")
print(f"  M1 vs M0 (Anillo solo):      Δ={nll_M0-nll_M1:+.3f}  log10(BF)={(nll_M0-nll_M1)/np.log(10):+.3f}")

# ============================================================
print(f"\n{'='*72}")
print(f"PARAMETROS MAP - MASAS CON FORMULA CORRECTA")
print(f"{'='*72}")

for name, p, nll_val in [('M2 (Kep solo)', p_M2, nll_M2),
                          ('M3 (Kep+ring)', p_M3, nll_M3)]:
    K_val = p[0]; P_val = p[1]; e_val = p[2]; om_val = p[3]
    A_r_val = p[5]; P_r_val = p[6]
    sigma_gp_val = np.exp(p[8])
    
    M_planet = msini_from_K_earth(K_val, P_val, M_star=0.972, e=e_val)
    
    print(f"\n  {name}:")
    print(f"    Kepleriano:")
    print(f"      K       = {K_val:.3f} m/s")
    print(f"      P       = {P_val:.3f} d")
    print(f"      e       = {e_val:.3f}")
    print(f"      omega   = {np.degrees(om_val):+.1f} deg")
    print(f"      M·sini  = {M_planet:.2f} M_Earth  (formula corregida)")
    if name == 'M3 (Kep+ring)':
        print(f"    Anillo:")
        print(f"      A       = {A_r_val:.3f} m/s")
        print(f"      P_ring  = {P_r_val:.2f} d")
    print(f"    GP actividad: sigma = {sigma_gp_val:.3f} m/s")

# Save
out = {
    'criterio_exclusion': f'std<{STD_THRESHOLD}, N>={MIN_NIGHTS}',
    'sub_inst_incluidos': sub_inst_names,
    'N_total': int(len(t)),
    'baseline_yr': float((t.max()-t.min())/365.25),
    'priors_v3_final': {
        'P_planet': '[217, 247] uniforme',
        'e_planet': '[0, 0.5] uniforme',
        'P_ring': '[400, 500] uniforme',
    },
    'M0_nll': float(nll_M0), 'M1_nll': float(nll_M1),
    'M2_nll': float(nll_M2), 'M3_nll': float(nll_M3),
    'M2_K': float(p_M2[0]), 'M2_P': float(p_M2[1]), 'M2_e': float(p_M2[2]),
    'M2_M_earth': float(msini_from_K_earth(p_M2[0], p_M2[1], 0.972, p_M2[2])),
    'M2_sigma_gp': float(np.exp(p_M2[8])),
    'M3_K': float(p_M3[0]), 'M3_P': float(p_M3[1]), 'M3_e': float(p_M3[2]),
    'M3_omega_deg': float(np.degrees(p_M3[3])),
    'M3_M_earth': float(msini_from_K_earth(p_M3[0], p_M3[1], 0.972, p_M3[2])),
    'M3_A_ring': float(p_M3[5]), 'M3_P_ring': float(p_M3[6]),
    'M3_sigma_gp': float(np.exp(p_M3[8])),
    'log10BF_M3_vs_M0': float((nll_M0-nll_M3)/np.log(10)),
    'log10BF_M2_vs_M0': float((nll_M0-nll_M2)/np.log(10)),
    'log10BF_M1_vs_M0': float((nll_M0-nll_M1)/np.log(10)),
    'log10BF_M3_vs_M2': float((nll_M2-nll_M3)/np.log(10)),
    'log10BF_M3_vs_M1': float((nll_M1-nll_M3)/np.log(10)),
}
with open(OUTPUT, 'w') as f:
    json.dump(out, f, indent=2)
print(f"\n[OK] {OUTPUT}")
