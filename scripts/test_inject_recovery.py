"""
Inject-recovery diagnostic test (V2 - fix walker independence).
Inyectamos Kepleriano sintetico (K=0.8, P=229.5, e=0.5) y vemos si MCMC lo recupera.
"""
import sys, os, json, time
sys.path.insert(0, '.')

import numpy as np
import pandas as pd
import emcee
import celerite2
from celerite2 import terms
from utils.rv_mass import msini_from_K_earth

INPUT  = "data_interim/qc/master_v3_v2_nightly_residuals_clean.csv"
MAP_JSON = "data_interim/qc/test2_map_v3_final.json"
OUTDIR = "data_interim/qc"

STD_THRESHOLD = 5.0
MIN_NIGHTS = 10
E_PRIOR_MAX = 0.85

K_INJ = 0.80
P_INJ = 229.5
E_INJ = 0.50
OMEGA_INJ = 0.0
T0_INJ = 2456000.0

def load_homogeneous_data():
    df = pd.read_csv(INPUT)
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    keep_subs = []
    for sub in sorted(df['sub_inst'].unique()):
        m = df['sub_inst'] == sub
        n = m.sum()
        std = df.loc[m, 'rv_resid_c'].std()
        if n >= MIN_NIGHTS and std <= STD_THRESHOLD:
            keep_subs.append(sub)
    df = df[df['sub_inst'].isin(keep_subs)].copy()
    df['rv_resid_c'] = df.groupby('sub_inst')['rv_resid'].transform(
        lambda x: x - x.median())
    df['inst_idx'] = df['sub_inst'].map(
        {s: i for i, s in enumerate(sorted(df['sub_inst'].unique()))})
    df = df.sort_values('bjd').reset_index(drop=True)
    return df, sorted(df['sub_inst'].unique())

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

def build_gp_kernel(sigma, period_rot):
    Q = 1.0 / np.sqrt(2)
    w0 = 2*np.pi/period_rot
    S0 = sigma**2 / (w0 * Q)
    return terms.SHOTerm(S0=S0, w0=w0, Q=Q)

# =====================================================================
# Modelo M2 puro: vector de 14 parametros + 2*n_inst (sin anillo)
# Layout: [K, P, e, omega, t0, log_sigma_gp, log_rho, log_prot, log_tau,
#          offset_0, ..., offset_{n_inst-1}, log_jitter_0, ..., log_jitter_{n_inst-1}]
# =====================================================================
def n_params_M2(n_inst):
    return 9 + 2*n_inst   # 5 Kep + 4 GP + 2*n_inst inst

def log_posterior_M2(params, t, rv, rv_err, inst_idx, n_inst):
    """M2 puro: solo Kep + GP + offsets/jitters. Sin anillo."""
    try:
        K = params[0]; P_p = params[1]; e = params[2]
        omega = params[3]; t0 = params[4]
        sigma_gp = np.exp(params[5])
        prot = np.exp(params[7])
        offsets = params[9:9+n_inst]
        jitters = np.exp(params[9+n_inst:9+2*n_inst])

        if K < 0 or K > 20: return -np.inf
        if P_p < 217 or P_p > 247: return -np.inf
        if e < 0 or e >= E_PRIOR_MAX: return -np.inf
        if abs(omega) > np.pi: return -np.inf
        if t0 < 2453000 or t0 > 2461000: return -np.inf
        if params[5] < -3 or params[5] > 4: return -np.inf
        if params[6] < 0 or params[6] > 7: return -np.inf
        lp = -0.5 * ((params[7] - np.log(38.65))/0.05)**2
        if params[7] < 3.0 or params[7] > 4.5: return -np.inf
        if params[8] < 0 or params[8] > 7: return -np.inf
        for off in offsets:
            if abs(off) > 30: return -np.inf
            lp += -0.5 * (off/15)**2
        for lj in params[9+n_inst:9+2*n_inst]:
            if lj < -3 or lj > 3: return -np.inf

        rv_model = keplerian_rv(t, K, P_p, e, omega, t0)
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

def init_walkers_M2(n_walkers, n_inst):
    n_dim = n_params_M2(n_inst)
    pos = np.zeros((n_walkers, n_dim))
    pos[:, 0]  = np.random.uniform(0.1, 5.0, n_walkers)         # K
    pos[:, 1]  = np.random.uniform(218, 246, n_walkers)         # P
    pos[:, 2]  = np.random.uniform(0.0, 0.80, n_walkers)        # e
    pos[:, 3]  = np.random.uniform(-np.pi, np.pi, n_walkers)    # omega
    pos[:, 4]  = np.random.uniform(2454000, 2459000, n_walkers) # t0
    pos[:, 5]  = np.random.uniform(-1, 2, n_walkers)            # log_sigma_gp
    pos[:, 6]  = np.random.uniform(2, 6, n_walkers)             # log_rho
    pos[:, 7]  = np.log(38.65) + 0.03*np.random.randn(n_walkers)# log_prot
    pos[:, 8]  = np.random.uniform(2, 6, n_walkers)             # log_tau
    for i in range(n_inst):
        pos[:, 9+i] = np.random.uniform(-5, 5, n_walkers)
        pos[:, 9+n_inst+i] = np.random.uniform(-2, 1, n_walkers)
    pos[:, 0] = np.clip(pos[:, 0], 1e-3, 19.9)
    pos[:, 2] = np.clip(pos[:, 2], 1e-3, E_PRIOR_MAX - 0.01)
    return pos

def run_mcmc_M2(t, rv, rv_err, inst_idx, n_inst, label, n_walkers=64, n_steps=20000):
    pos = init_walkers_M2(n_walkers, n_inst)
    n_dim = n_params_M2(n_inst)
    
    sampler = emcee.EnsembleSampler(
        n_walkers, n_dim, log_posterior_M2,
        args=(t, rv, rv_err, inst_idx, n_inst),
    )
    
    print(f"\n[{label}] MCMC: {n_walkers}w x {n_steps}st  n_dim={n_dim}")
    t0 = time.time()
    for sample in sampler.sample(pos, iterations=n_steps, progress=False):
        i = sampler.iteration
        if i % 2000 == 0:
            elapsed = time.time() - t0
            try:
                tau = sampler.get_autocorr_time(tol=0)
                tau_mean = np.mean(tau)
                print(f"  [{label}] step {i}/{n_steps} t={elapsed/60:.1f}m tau={tau_mean:.0f}")
            except:
                print(f"  [{label}] step {i}/{n_steps} t={elapsed/60:.1f}m")
    
    elapsed_total = time.time() - t0
    print(f"  [{label}] DONE in {elapsed_total/60:.1f}m")
    
    try:
        tau_final = sampler.get_autocorr_time(tol=0)
        burnin = int(3 * np.max(tau_final))
        thin = max(int(0.5 * np.min(tau_final)), 1)
    except Exception:
        burnin = sampler.iteration // 4
        thin = 1
    
    chain = sampler.get_chain(discard=burnin, thin=thin, flat=True)
    log_prob = sampler.get_log_prob(discard=burnin, thin=thin, flat=True)
    return chain, log_prob

# =====================================================================
def main():
    print(f"="*72)
    print(f"INJECT-RECOVERY DIAGNOSTIC TEST (V2 - walker independence fix)")
    print(f"  Injection: K={K_INJ}, P={P_INJ}, e={E_INJ}, omega={OMEGA_INJ}")
    print(f"="*72)
    
    df, sub_inst_names = load_homogeneous_data()
    n_inst = len(sub_inst_names)
    t = df['bjd'].values
    rv_orig = df['rv_resid_c'].values
    rv_err = np.maximum(df['rv_err'].values, 0.1)
    inst_idx = df['inst_idx'].values
    
    print(f"\n  Data: N={len(t)}, n_inst={n_inst}")
    
    # CONTROL
    print(f"\n{'='*72}")
    print(f"CONTROL: MCMC sobre datos ORIGINALES (sin inyeccion)")
    print(f"{'='*72}")
    chain_ctrl, lp_ctrl = run_mcmc_M2(t, rv_orig, rv_err, inst_idx, n_inst, "CONTROL")
    
    K_ctrl = np.median(chain_ctrl[:, 0])
    K_ctrl_lo, K_ctrl_hi = np.percentile(chain_ctrl[:, 0], [16, 84])
    K_ctrl_sigma = K_ctrl / ((K_ctrl_hi - K_ctrl_lo) / 2)
    e_ctrl = np.median(chain_ctrl[:, 2])
    P_ctrl = np.median(chain_ctrl[:, 1])
    
    print(f"\n  CONTROL posteriors:")
    print(f"    K = {K_ctrl:.3f} (+{K_ctrl_hi-K_ctrl:.3f} / -{K_ctrl-K_ctrl_lo:.3f}) m/s")
    print(f"    P = {P_ctrl:.3f} d")
    print(f"    e = {e_ctrl:.3f}")
    print(f"    K_sigma = {K_ctrl_sigma:.2f}")
    
    # INJECTION
    print(f"\n{'='*72}")
    print(f"INJECTION: a\u00f1adir Kepleriano K={K_INJ}, P={P_INJ}, e={E_INJ}")
    print(f"{'='*72}")
    
    rv_injected_signal = keplerian_rv(t, K_INJ, P_INJ, E_INJ, OMEGA_INJ, T0_INJ)
    rv_inj = rv_orig + rv_injected_signal
    
    print(f"  Senal inyectada: amplitud peak-to-peak = {rv_injected_signal.max()-rv_injected_signal.min():.3f} m/s")
    print(f"  Std de la inyeccion: {rv_injected_signal.std():.3f} m/s")
    
    chain_inj, lp_inj = run_mcmc_M2(t, rv_inj, rv_err, inst_idx, n_inst, "INJECTED")
    
    K_inj_rec = np.median(chain_inj[:, 0])
    K_inj_lo, K_inj_hi = np.percentile(chain_inj[:, 0], [16, 84])
    K_inj_sigma = K_inj_rec / ((K_inj_hi - K_inj_lo) / 2)
    e_inj_rec = np.median(chain_inj[:, 2])
    P_inj_rec = np.median(chain_inj[:, 1])
    
    print(f"\n  INJECTED posteriors (esperado: K\u22480.8, P\u2248229.5, e\u22480.50):")
    print(f"    K = {K_inj_rec:.3f} (+{K_inj_hi-K_inj_rec:.3f} / -{K_inj_rec-K_inj_lo:.3f}) m/s")
    print(f"    P = {P_inj_rec:.3f} d")
    print(f"    e = {e_inj_rec:.3f}")
    print(f"    K_sigma = {K_inj_sigma:.2f}")
    
    # DIAGNOSTICO
    print(f"\n{'='*72}")
    print(f"DIAGNOSTICO INJECT-RECOVERY")
    print(f"{'='*72}")
    
    K_diff = K_inj_rec - K_ctrl
    
    print(f"\n  Comparacion K:")
    print(f"    K_INJ (inyectado):              {K_INJ:.3f}")
    print(f"    K_ctrl (sin inyectar):          {K_ctrl:.3f}  K_sigma={K_ctrl_sigma:.2f}")
    print(f"    K_inj_rec (con inyeccion):      {K_inj_rec:.3f}  K_sigma={K_inj_sigma:.2f}")
    print(f"    Diferencia (recovery - ctrl):   {K_diff:+.3f}  (esperado ~+{K_INJ:.2f})")
    
    print(f"\n  Comparacion e:")
    print(f"    e_INJ (inyectado):              {E_INJ:.3f}")
    print(f"    e_ctrl (sin inyectar):          {e_ctrl:.3f}")
    print(f"    e_inj_rec (con):                {e_inj_rec:.3f}")
    
    print(f"\n  VEREDICTO:")
    if K_inj_sigma > 3 and abs(K_diff - K_INJ) < 0.3:
        print(f"  [SE\u00d1AL DETECTABLE] Dataset recupera K={K_INJ} con K_sigma={K_inj_sigma:.1f}")
        print(f"    Que el control de {K_ctrl:.2f} m/s con K_sigma={K_ctrl_sigma:.1f} es MARGINAL,")
        print(f"    sugiere que la senal control es debil pero podria ser real.")
    elif K_inj_sigma > 2:
        print(f"  [DETECCION DEBIL] Recuperacion parcial K_sigma={K_inj_sigma:.1f}")
        print(f"    Dataset cerca del limite de detectabilidad.")
    else:
        print(f"  [LIMITE FUNDAMENTAL] No se recupera la inyeccion (K_sigma={K_inj_sigma:.1f})")
        print(f"    El dataset NO puede detectar K={K_INJ} con estos parametros.")
        print(f"    Cualquier 'senal' del control es indistinguible de ruido.")
        print(f"    SERVAL ES ESENCIAL para conclusion mas alla.")
    
    # Save
    np.save(f"{OUTDIR}/inject_recovery_control_chain.npy", chain_ctrl)
    np.save(f"{OUTDIR}/inject_recovery_injected_chain.npy", chain_inj)
    
    summary = {
        'injection_params': {
            'K_inj': K_INJ, 'P_inj': P_INJ, 'e_inj': E_INJ,
            'omega_inj': OMEGA_INJ, 't0_inj': T0_INJ,
        },
        'control': {
            'K_median': float(K_ctrl), 'K_16': float(K_ctrl_lo), 'K_84': float(K_ctrl_hi),
            'P_median': float(P_ctrl), 'e_median': float(e_ctrl),
            'K_sigma': float(K_ctrl_sigma),
            'log_prob_max': float(np.max(lp_ctrl)),
        },
        'injected': {
            'K_median': float(K_inj_rec), 'K_16': float(K_inj_lo), 'K_84': float(K_inj_hi),
            'P_median': float(P_inj_rec), 'e_median': float(e_inj_rec),
            'K_sigma': float(K_inj_sigma),
            'log_prob_max': float(np.max(lp_inj)),
        },
        'diagnostic': {
            'K_diff_recovered_minus_control': float(K_diff),
            'K_diff_expected': float(K_INJ),
            'K_recovery_fraction': float(K_diff / K_INJ),
        },
    }
    with open(f"{OUTDIR}/inject_recovery_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n[OK] {OUTDIR}/inject_recovery_summary.json")

if __name__ == "__main__":
    main()
