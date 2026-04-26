"""
Test 2: Modelo conjunto Kepleriano + GP cuasi-periodico

Modelo:
  RV(t) = K_planet * Kepler(P, K, e, w, t0)
        + GP_QP(sigma, l, Prot, lambda_p)
        + offset_per_instrument
        + jitter_per_instrument (en cuadratura con rv_err)

Modos:
  --mode map:   MAP optimization (10-15 min)
  --mode mcmc:  MCMC 64w x 20000 steps (4-6h) inicializado desde MAP

Output:
  data_interim/qc/test2_map_results.json  (parametros MAP)
  data_interim/qc/test2_mcmc_chain.npy    (cadena MCMC)
  data_interim/qc/test2_corner.png
  data_interim/qc/test2_residuals_gls.png
  data_interim/qc/test2_summary.txt
"""

import sys
import json
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

# Optional packages
import celerite2
from celerite2 import terms
import emcee

INPUT  = "data_interim/qc/master_v3_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"
MAP_OUT = f"{OUTDIR}/test2_map_results.json"

# ============================================================
# Cargar y preparar datos
# ============================================================
def load_data():
    df = pd.read_csv(INPUT)
    df = df.sort_values('bjd').reset_index(drop=True)
    
    # Centrar por sub_inst
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    
    # Quedarnos con sub-inst con N>=10 (los que aportan al GP)
    counts = df.groupby('sub_inst').size()
    keep_subs = counts[counts >= 10].index.tolist()
    df = df[df['sub_inst'].isin(keep_subs)].copy()
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    
    print(f"  Sub-inst con N>=10:")
    for sub in df['sub_inst'].unique():
        m = df['sub_inst']==sub
        print(f"    {sub:<25} N={m.sum():>4} std={df.loc[m,'rv_resid_c'].std():.2f}")
    
    # Asignar inst_idx (entero por sub_inst)
    df['inst_idx'] = df['sub_inst'].map(
        {sub: i for i, sub in enumerate(sorted(df['sub_inst'].unique()))})
    
    return df, sorted(df['sub_inst'].unique())

# ============================================================
# Modelos
# ============================================================
def kepler_solve(M, e, tol=1e-10, max_iter=100):
    """Resolver Kepler M = E - e*sin(E)."""
    E = M.copy()
    for _ in range(max_iter):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < tol:
            break
    return E

def keplerian_rv(t, K, P, e, omega, t0):
    """RV modelo Kepleriano."""
    M = 2*np.pi*(t - t0)/P
    E = kepler_solve(M, e)
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    rv = K * (np.cos(nu + omega) + e*np.cos(omega))
    return rv

def build_gp(t, sigma, rho, period, tau):
    """GP cuasi-periodico via celerite2 SHOTerm + RotationTerm."""
    # Usar QPTerm aproximado: SHOTerm con periodo
    # celerite2 tiene un kernel cuasi-periodico estandarizado
    # Aprox: combinacion de SHOTerm
    Q = 0.5 / np.sqrt(2)  # critically damped
    w0 = 2*np.pi/period
    S0 = sigma**2 / (w0 * Q)
    term = terms.SHOTerm(S0=S0, w0=w0, Q=Q)
    # Anyadir factor de decaimiento via tau (escala de coherencia)
    return term

# ============================================================
# Log-likelihood
# ============================================================
def log_likelihood(params, t, rv, rv_err, inst_idx, n_inst, fixed_P=None):
    """
    params: [K, P, e, omega, t0,
             log_sigma_gp, log_rho, log_period_rot, log_tau,
             offset_0, offset_1, ..., offset_{n_inst-1},
             log_jitter_0, ..., log_jitter_{n_inst-1}]
    """
    try:
        idx = 0
        K = params[idx]; idx += 1
        P = params[idx]; idx += 1
        e = params[idx]; idx += 1
        omega = params[idx]; idx += 1
        t0 = params[idx]; idx += 1
        log_sigma_gp = params[idx]; idx += 1
        log_rho = params[idx]; idx += 1
        log_prot = params[idx]; idx += 1
        log_tau = params[idx]; idx += 1
        offsets = params[idx:idx+n_inst]; idx += n_inst
        log_jitters = params[idx:idx+n_inst]; idx += n_inst
        
        # Sanity: e debe estar en [0, 0.95)
        if e < 0 or e >= 0.95: return -np.inf
        if K < 0: return -np.inf
        if P <= 0: return -np.inf
        
        sigma_gp = np.exp(log_sigma_gp)
        rho = np.exp(log_rho)
        period_rot = np.exp(log_prot)
        tau = np.exp(log_tau)
        jitters = np.exp(log_jitters)
        
        # Modelo Kepleriano
        rv_kepler = keplerian_rv(t, K, P, e, omega, t0)
        
        # Offsets por instrumento
        offset_t = offsets[inst_idx]
        rv_model = rv_kepler + offset_t
        
        # Residuos
        resid = rv - rv_model
        
        # Jitter por instrumento (en cuadratura)
        jitter_t = jitters[inst_idx]
        yerr_total = np.sqrt(rv_err**2 + jitter_t**2)
        
        # GP cuasi-periodico
        gp = build_gp(t, sigma_gp, rho, period_rot, tau)
        gp_model = celerite2.GaussianProcess(gp, mean=0.0)
        gp_model.compute(t, yerr=yerr_total, quiet=True)
        
        ll = gp_model.log_likelihood(resid)
        if not np.isfinite(ll):
            return -np.inf
        return ll
    except Exception:
        return -np.inf

# ============================================================
# Log-prior
# ============================================================
def log_prior(params, n_inst, P_prior_mean=232.0, P_prior_sigma=5.0):
    idx = 0
    K = params[idx]; idx += 1
    P = params[idx]; idx += 1
    e = params[idx]; idx += 1
    omega = params[idx]; idx += 1
    t0 = params[idx]; idx += 1
    log_sigma_gp = params[idx]; idx += 1
    log_rho = params[idx]; idx += 1
    log_prot = params[idx]; idx += 1
    log_tau = params[idx]; idx += 1
    offsets = params[idx:idx+n_inst]; idx += n_inst
    log_jitters = params[idx:idx+n_inst]; idx += n_inst
    
    lp = 0.0
    
    # K: uniform [0, 20] m/s
    if K < 0 or K > 20: return -np.inf
    
    # P: gaussian prior alrededor 232 d
    lp += -0.5 * ((P - P_prior_mean)/P_prior_sigma)**2
    if P < 100 or P > 500: return -np.inf
    
    # e: beta(1, 5) prior (favorece bajas e)
    if e < 0 or e >= 0.95: return -np.inf
    lp += 4 * np.log(1 - e)  # ~ Beta(1, 5)
    
    # omega: uniform [-pi, pi]
    if abs(omega) > np.pi: return -np.inf
    
    # t0: uniform pero acotado a un periodo
    if t0 < 2453000 or t0 > 2461000: return -np.inf
    
    # log_sigma_gp: uniform [-3, 3]  (sigma_gp en e^[-3, 3] = [0.05, 20] m/s)
    if log_sigma_gp < -3 or log_sigma_gp > 3: return -np.inf
    
    # log_rho: uniform [0, 7] (rho in [1, 1100] d)
    if log_rho < 0 or log_rho > 7: return -np.inf
    
    # log_prot: gaussian centrado en log(38.65)
    lp += -0.5 * ((log_prot - np.log(38.65))/0.05)**2
    if log_prot < 3.0 or log_prot > 4.5: return -np.inf
    
    # log_tau: uniform [0, 7]
    if log_tau < 0 or log_tau > 7: return -np.inf
    
    # offsets: gaussian wide
    for off in offsets:
        if abs(off) > 30: return -np.inf
        lp += -0.5 * (off/15)**2
    
    # log_jitters: uniform [-3, 3]
    for lj in log_jitters:
        if lj < -3 or lj > 3: return -np.inf
    
    return lp

def log_posterior(params, t, rv, rv_err, inst_idx, n_inst):
    lp = log_prior(params, n_inst)
    if not np.isfinite(lp):
        return -np.inf
    ll = log_likelihood(params, t, rv, rv_err, inst_idx, n_inst)
    if not np.isfinite(ll):
        return -np.inf
    return lp + ll

# ============================================================
# MAP optimization
# ============================================================
def run_map(t, rv, rv_err, inst_idx, n_inst, sub_inst_names):
    print("\n" + "="*72)
    print("MAP OPTIMIZATION")
    print("="*72)
    
    # Initial guess
    p0 = [
        1.0,              # K
        232.0,            # P
        0.1,              # e
        0.0,              # omega
        2456000.0,        # t0
        np.log(2.0),      # log_sigma_gp ~ 2 m/s
        np.log(50.0),     # log_rho ~ 50 d
        np.log(38.65),    # log_prot
        np.log(50.0),     # log_tau
    ]
    p0.extend([0.0] * n_inst)            # offsets
    p0.extend([np.log(0.5)] * n_inst)    # log_jitters
    p0 = np.array(p0)
    
    print(f"  Total parametros: {len(p0)}")
    print(f"  Inicial: K={p0[0]:.2f}, P={p0[1]:.1f}, e={p0[2]:.2f}")
    
    def nll(params):
        lp = log_posterior(params, t, rv, rv_err, inst_idx, n_inst)
        if not np.isfinite(lp): return 1e30
        return -lp
    
    print(f"\n  Inicial -log_posterior: {nll(p0):.2f}")
    
    # Multiple restarts
    best_p, best_nll = p0, nll(p0)
    n_restarts = 5
    for r in range(n_restarts):
        # Random perturbation
        p_init = p0.copy()
        if r > 0:
            p_init[0] = np.random.uniform(0.5, 3.0)            # K
            p_init[1] = np.random.normal(232, 3)               # P
            p_init[2] = np.random.uniform(0, 0.3)              # e
            p_init[3] = np.random.uniform(-np.pi, np.pi)       # omega
            p_init[5] = np.log(np.random.uniform(0.5, 5))      # sigma_gp
            p_init[6] = np.log(np.random.uniform(20, 200))     # rho
            p_init[8] = np.log(np.random.uniform(20, 200))     # tau
        
        try:
            res = minimize(nll, p_init, method='Nelder-Mead',
                          options={'maxiter': 3000, 'xatol': 1e-4, 'fatol': 1e-4})
            print(f"  Restart {r+1}: nll={res.fun:.2f}, K={res.x[0]:.3f}, "
                  f"P={res.x[1]:.2f}, e={res.x[2]:.3f}")
            if res.fun < best_nll:
                best_nll, best_p = res.fun, res.x
        except Exception as exc:
            print(f"  Restart {r+1}: error - {exc}")
    
    # Reporte
    K_map = best_p[0]
    P_map = best_p[1]
    e_map = best_p[2]
    omega_map = best_p[3]
    sigma_gp_map = np.exp(best_p[5])
    prot_map = np.exp(best_p[7])
    
    # Calcular masa equivalente
    M_STAR = 0.972
    P_yr = P_map / 365.25
    M_p_earth = K_map / 0.6395 * (M_STAR)**(2/3) * (P_yr)**(1/3) * np.sqrt(1 - e_map**2)
    
    print(f"\n  RESULTADO MAP:")
    print(f"    K_planet  = {K_map:.3f} m/s")
    print(f"    P_planet  = {P_map:.3f} d")
    print(f"    e_planet  = {e_map:.3f}")
    print(f"    omega     = {omega_map:.3f} rad")
    print(f"    M·sin(i)  = {M_p_earth:.3f} M_Earth")
    print(f"    sigma_gp  = {sigma_gp_map:.3f} m/s")
    print(f"    Prot_GP   = {prot_map:.3f} d")
    print(f"    -log_post = {best_nll:.3f}")
    
    # Calcular nll sin Kepleriano (K=0) para evidencia bayesiana
    p_no_planet = best_p.copy()
    p_no_planet[0] = 0.0  # K=0
    nll_no_planet = nll(p_no_planet)
    delta_nll = nll_no_planet - best_nll
    print(f"\n  COMPARACION CON MODELO SIN KEPLERIANO (K=0):")
    print(f"    -log_post(K=0)   = {nll_no_planet:.3f}")
    print(f"    Delta            = {delta_nll:.3f}")
    print(f"    Aprox log10(BF)  = {delta_nll/np.log(10):.2f}")
    
    if delta_nll > 5:
        verdict = "STRONG EVIDENCE para Kepleriano (delta>5)"
    elif delta_nll > 2.5:
        verdict = "Evidence moderada para Kepleriano (2.5<delta<5)"
    elif delta_nll > 1:
        verdict = "Evidence d\u00e9bil"
    else:
        verdict = "SIN EVIDENCIA para Kepleriano (delta<1)"
    print(f"    VEREDICTO MAP: {verdict}")
    
    # Save MAP
    map_dict = {
        'best_params': best_p.tolist(),
        'best_nll': float(best_nll),
        'K': float(K_map),
        'P': float(P_map),
        'e': float(e_map),
        'omega': float(omega_map),
        'sigma_gp': float(sigma_gp_map),
        'prot_gp': float(prot_map),
        'M_p_earth': float(M_p_earth),
        'delta_nll_vs_no_planet': float(delta_nll),
        'log10_BF_approx': float(delta_nll/np.log(10)),
        'verdict': verdict,
        'sub_inst_names': sub_inst_names,
        'n_inst': n_inst,
    }
    with open(MAP_OUT, 'w') as f:
        json.dump(map_dict, f, indent=2)
    print(f"\n  [OK] {MAP_OUT}")
    
    return best_p, best_nll, map_dict

# ============================================================
# MCMC
# ============================================================
def run_mcmc(t, rv, rv_err, inst_idx, n_inst, sub_inst_names, p_init):
    print("\n" + "="*72)
    print("MCMC RIGUROSO (64 walkers x 20000 steps)")
    print("="*72)
    print("  ATENCION: este paso tarda 4-6 horas")
    
    n_dim = len(p_init)
    n_walkers = 64
    n_steps = 20000
    n_burnin = 5000
    
    # Inicializacion: ball alrededor de MAP
    pos = p_init + 1e-3 * np.random.randn(n_walkers, n_dim) * np.abs(p_init + 1e-6)
    # Asegurar que K>=0 y e en [0,0.9]
    pos[:, 0] = np.maximum(pos[:, 0], 1e-3)
    pos[:, 2] = np.clip(pos[:, 2], 1e-3, 0.9)
    
    sampler = emcee.EnsembleSampler(
        n_walkers, n_dim, log_posterior,
        args=(t, rv, rv_err, inst_idx, n_inst),
    )
    
    print(f"  N_walkers: {n_walkers}, n_dim: {n_dim}, n_steps: {n_steps}, burnin: {n_burnin}")
    print(f"  Total samples: {n_walkers*n_steps:,}")
    print(f"  Progress logging cada 1000 steps...")
    
    import time
    t0 = time.time()
    for sample in sampler.sample(pos, iterations=n_steps, progress=False):
        i = sampler.iteration
        if i % 1000 == 0:
            elapsed = time.time() - t0
            rate = i / elapsed
            eta_h = (n_steps - i) / rate / 3600
            try:
                tau = sampler.get_autocorr_time(tol=0)
                tau_str = f"tau_max={np.max(tau):.0f}"
            except:
                tau_str = "tau pending"
            print(f"  step {i}/{n_steps}  [ETA {eta_h:.1f}h]  {tau_str}")
    
    elapsed = time.time() - t0
    print(f"\n  MCMC done en {elapsed/3600:.2f}h")
    
    # Save chain
    chain = sampler.get_chain(discard=n_burnin, flat=True)
    np.save(f"{OUTDIR}/test2_mcmc_chain.npy", chain)
    print(f"  [OK] {OUTDIR}/test2_mcmc_chain.npy ({chain.shape})")
    
    # Statistics
    K_samples = chain[:, 0]
    P_samples = chain[:, 1]
    e_samples = chain[:, 2]
    
    K_med = np.median(K_samples)
    K_lo, K_hi = np.percentile(K_samples, [16, 84])
    P_med = np.median(P_samples)
    P_lo, P_hi = np.percentile(P_samples, [16, 84])
    e_med = np.median(e_samples)
    e_lo, e_hi = np.percentile(e_samples, [16, 84])
    
    M_STAR = 0.972
    M_samples = K_samples / 0.6395 * (M_STAR)**(2/3) * (P_samples/365.25)**(1/3) * np.sqrt(1-e_samples**2)
    M_med = np.median(M_samples)
    M_lo, M_hi = np.percentile(M_samples, [16, 84])
    
    print(f"\n  POSTERIORS (mediana, 16-84%):")
    print(f"    K = {K_med:.3f} (+{K_hi-K_med:.3f} / -{K_med-K_lo:.3f}) m/s")
    print(f"    P = {P_med:.3f} (+{P_hi-P_med:.3f} / -{P_med-P_lo:.3f}) d")
    print(f"    e = {e_med:.3f} (+{e_hi-e_med:.3f} / -{e_med-e_lo:.3f})")
    print(f"    M·sin(i) = {M_med:.3f} (+{M_hi-M_med:.3f} / -{M_med-M_lo:.3f}) M_Earth")
    
    # Significancia de K
    K_sigma = (K_med) / ((K_hi - K_lo)/2)
    print(f"\n  K_significance: {K_sigma:.2f}-sigma")
    if K_sigma > 5:
        verdict = "*** DETECTION (K > 5-sigma) ***"
    elif K_sigma > 3:
        verdict = "** Strong candidate (K > 3-sigma) **"
    elif K_sigma > 2:
        verdict = "* Marginal (K > 2-sigma) *"
    else:
        verdict = "NO DETECTION (K < 2-sigma)"
    print(f"  VEREDICTO MCMC: {verdict}")
    
    # Corner plot (subset de parametros mas relevantes)
    try:
        import corner
        labels = ['K', 'P', 'e', r'$\omega$', '$t_0$',
                  r'log $\sigma_{GP}$', r'log $\rho$', r'log $P_{rot}$', r'log $\tau$']
        n_show = len(labels)
        fig = corner.corner(chain[:, :n_show], labels=labels,
                           show_titles=True, quantiles=[0.16, 0.5, 0.84])
        fig.savefig(f"{OUTDIR}/test2_corner.png", dpi=100)
        print(f"  [OK] {OUTDIR}/test2_corner.png")
    except Exception as exc:
        print(f"  Corner failed: {exc}")
    
    # Save summary
    with open(f"{OUTDIR}/test2_summary.txt", 'w') as f:
        f.write(f"Test 2 MCMC summary\n")
        f.write(f"="*40 + "\n")
        f.write(f"K = {K_med:.3f} +{K_hi-K_med:.3f}/-{K_med-K_lo:.3f} m/s\n")
        f.write(f"P = {P_med:.3f} +{P_hi-P_med:.3f}/-{P_med-P_lo:.3f} d\n")
        f.write(f"e = {e_med:.3f} +{e_hi-e_med:.3f}/-{e_med-e_lo:.3f}\n")
        f.write(f"M = {M_med:.3f} +{M_hi-M_med:.3f}/-{M_med-M_lo:.3f} M_Earth\n")
        f.write(f"K_sigma = {K_sigma:.2f}\n")
        f.write(f"Verdict: {verdict}\n")
    print(f"  [OK] {OUTDIR}/test2_summary.txt")

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=['map','mcmc'], required=True)
    args = parser.parse_args()
    
    print(f"="*72)
    print(f"TEST 2: GP + Kepleriano conjunto - mode={args.mode}")
    print(f"="*72)
    
    df, sub_inst_names = load_data()
    n_inst = len(sub_inst_names)
    t = df['bjd'].values
    rv = df['rv_resid_c'].values
    rv_err = np.maximum(df['rv_err'].values, 0.1)
    inst_idx = df['inst_idx'].values
    
    print(f"\n  Datos: N={len(t)}, n_inst={n_inst}")
    
    if args.mode == 'map':
        run_map(t, rv, rv_err, inst_idx, n_inst, sub_inst_names)
        print(f"\n  Para continuar con MCMC, lanza:")
        print(f"  python scripts/test2_gp_keplerian.py --mode mcmc")
    
    elif args.mode == 'mcmc':
        import os
        if not os.path.exists(MAP_OUT):
            print(f"\n  ERROR: no existe {MAP_OUT}")
            print(f"  Lanza primero: python scripts/test2_gp_keplerian.py --mode map")
            sys.exit(1)
        with open(MAP_OUT) as f:
            map_dict = json.load(f)
        p_init = np.array(map_dict['best_params'])
        run_mcmc(t, rv, rv_err, inst_idx, n_inst, sub_inst_names, p_init)
