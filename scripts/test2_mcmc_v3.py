"""
Test 2 V3 MCMC: 4 modelos (M0/M1/M2/M3) con priors V3 + e abierto a 0.7

Modelos:
  M0: GP + offsets (nulo)
  M1: GP + ring (anillo solo)
  M2: GP + Kepleriano (planet solo) <- LOG10BF=+1.09 en MAP V3
  M3: GP + Kepleriano + ring (completo)

Configuracion:
  64 walkers x 20000 steps (production)
  Convergencia tau-based: para cuando tau cambia <1% en 100 steps
  Inicializacion: --init {map5, prior}
    map5 = MAP +/- 5% (default, convergencia rapida ~1-2h por modelo)
    prior = uniforme en prior amplio (~3-6h por modelo)

Priors:
  P_planet uniforme [217, 247]
  e_planet uniforme [0, 0.7]   <- AMPLIO, permite chequear si MAP corta
  omega uniforme [-pi, pi]
  P_ring uniforme [400, 500]

Output:
  data_interim/qc/test2_mcmc_v3_M{0,1,2,3}.npy   (cadenas)
  data_interim/qc/test2_mcmc_v3_summary.json
  data_interim/qc/test2_mcmc_v3_corner.png       (M2 y M3)
  data_interim/qc/test2_mcmc_v3_progress.txt     (log de progreso)
"""
import sys, os, json, time, argparse
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

# ============================================================
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
        sigma_gp = np.exp(params[8])
        prot = np.exp(params[10])
        offsets = params[12:12+n_inst]
        jitters = np.exp(params[12+n_inst:12+2*n_inst])

        # PRIORS V3 + e amplio a 0.7
        if K < 0 or K > 20: return -np.inf
        if P_p < 217 or P_p > 247: return -np.inf
        if e < 0 or e >= 0.7: return -np.inf       # <- ABIERTO
        if abs(omega) > np.pi: return -np.inf
        if t0 < 2453000 or t0 > 2461000: return -np.inf
        if A_ring < 0 or A_ring > 10: return -np.inf
        if P_ring < 400 or P_ring > 500: return -np.inf
        if abs(phi_ring) > np.pi: return -np.inf
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
def init_walkers(n_walkers, n_dim, init_mode, p_map=None):
    """Inicializar walkers segun modo."""
    if init_mode == 'map5' and p_map is not None:
        # MAP +/- 5%
        scale = np.abs(p_map) * 0.05 + 1e-4
        pos = p_map[None,:] + scale[None,:] * np.random.randn(n_walkers, n_dim)
    else:
        # Prior amplio
        pos = np.zeros((n_walkers, n_dim))
        pos[:, 0]  = np.random.uniform(0.1, 5.0, n_walkers)        # K
        pos[:, 1]  = np.random.uniform(218, 246, n_walkers)        # P_planet
        pos[:, 2]  = np.random.uniform(0.0, 0.65, n_walkers)       # e
        pos[:, 3]  = np.random.uniform(-np.pi, np.pi, n_walkers)   # omega
        pos[:, 4]  = np.random.uniform(2454000, 2459000, n_walkers)# t0
        pos[:, 5]  = np.random.uniform(0.1, 5.0, n_walkers)        # A_ring
        pos[:, 6]  = np.random.uniform(405, 495, n_walkers)        # P_ring
        pos[:, 7]  = np.random.uniform(-np.pi, np.pi, n_walkers)   # phi_ring
        pos[:, 8]  = np.random.uniform(-1, 2, n_walkers)           # log_sigma_gp
        pos[:, 9]  = np.random.uniform(2, 6, n_walkers)            # log_rho
        pos[:, 10] = np.log(38.65) + 0.03*np.random.randn(n_walkers)
        pos[:, 11] = np.random.uniform(2, 6, n_walkers)            # log_tau
        # offsets
        n_inst = (n_dim - 12) // 2
        for i in range(n_inst):
            pos[:, 12+i] = np.random.uniform(-5, 5, n_walkers)
            pos[:, 12+n_inst+i] = np.random.uniform(-2, 1, n_walkers)
    # Asegurar rangos validos para K, e
    pos[:, 0] = np.clip(pos[:, 0], 1e-3, 19.9)
    pos[:, 2] = np.clip(pos[:, 2], 1e-3, 0.69)
    pos[:, 5] = np.clip(pos[:, 5], 1e-3, 9.9)
    return pos

# ============================================================
def get_map_for_model(map_data, model_name, n_inst):
    """Construir vector de parametros MAP para un modelo dado."""
    n_p = n_params(n_inst)
    p = np.zeros(n_p)
    if model_name == 'M2':
        p[0] = map_data['M2_K']
        p[1] = map_data['M2_P']
        p[2] = map_data['M2_e']
        p[3] = 0.0
        p[8] = np.log(map_data['M2_sigma_gp'])
    elif model_name == 'M3':
        p[0] = map_data['M3_K']
        p[1] = map_data['M3_P']
        p[2] = map_data['M3_e']
        p[3] = np.radians(map_data['M3_omega_deg'])
        p[5] = map_data['M3_A_ring']
        p[6] = map_data['M3_P_ring']
        p[8] = np.log(map_data['M3_sigma_gp'])
    else:
        # M0, M1: usar valores neutros
        p[0] = 1.0; p[1] = 232.0; p[2] = 0.2
        p[5] = 1.0; p[6] = 450.0
        p[8] = np.log(2.0)
    p[4] = 2456000.0
    p[7] = 0.0
    p[9] = np.log(50.0)
    p[10] = np.log(38.65)
    p[11] = np.log(50.0)
    p[12:12+n_inst] = 0.0
    p[12+n_inst:] = np.log(0.5)
    return p

# ============================================================
def run_mcmc_model(t, rv, rv_err, inst_idx, n_inst,
                   use_planet, use_ring, model_name,
                   n_walkers, n_steps, init_mode, p_map=None,
                   logfile=None):
    n_dim = n_params(n_inst)
    
    pos = init_walkers(n_walkers, n_dim, init_mode, p_map)
    
    sampler = emcee.EnsembleSampler(
        n_walkers, n_dim, log_posterior,
        args=(t, rv, rv_err, inst_idx, n_inst),
        kwargs={'use_planet': use_planet, 'use_ring': use_ring},
    )
    
    msg = (f"\n[{model_name}] MCMC: {n_walkers}w x {n_steps}st  "
           f"init={init_mode}  planet={use_planet}  ring={use_ring}")
    print(msg)
    if logfile: logfile.write(msg + "\n"); logfile.flush()
    
    t0 = time.time()
    autocorr_history = []
    converged_at = None
    last_tau = None
    
    for sample in sampler.sample(pos, iterations=n_steps, progress=False):
        i = sampler.iteration
        if i % 500 == 0:
            elapsed = time.time() - t0
            try:
                tau = sampler.get_autocorr_time(tol=0)
                tau_mean = np.mean(tau)
                autocorr_history.append(tau_mean)
                # Convergencia: tau cambia <1% en ultimos 5 chequeos
                if len(autocorr_history) > 5:
                    recent = autocorr_history[-5:]
                    if (np.std(recent) / np.mean(recent)) < 0.01 and i >= 50*tau_mean:
                        converged_at = i
                        msg = f"  [{model_name}] step {i}: CONVERGED tau={tau_mean:.0f}"
                        print(msg)
                        if logfile: logfile.write(msg + "\n"); logfile.flush()
                        break
                msg = f"  [{model_name}] step {i}/{n_steps} t={elapsed/60:.1f}m tau={tau_mean:.0f}"
            except Exception:
                msg = f"  [{model_name}] step {i}/{n_steps} t={elapsed/60:.1f}m tau=pending"
            print(msg)
            if logfile: logfile.write(msg + "\n"); logfile.flush()
    
    elapsed_total = time.time() - t0
    final_iter = sampler.iteration
    msg = f"  [{model_name}] DONE in {elapsed_total/60:.1f}m at step {final_iter}"
    print(msg)
    if logfile: logfile.write(msg + "\n"); logfile.flush()
    
    # Determinar burnin: 3*tau si converged, 25% si no
    try:
        tau_final = sampler.get_autocorr_time(tol=0)
        burnin = int(3 * np.max(tau_final))
        thin = int(0.5 * np.min(tau_final))
        thin = max(thin, 1)
    except Exception:
        burnin = final_iter // 4
        thin = 1
    
    chain = sampler.get_chain(discard=burnin, thin=thin, flat=True)
    log_prob = sampler.get_log_prob(discard=burnin, thin=thin, flat=True)
    
    msg = (f"  [{model_name}] burnin={burnin} thin={thin} "
           f"final_chain={chain.shape}")
    print(msg)
    if logfile: logfile.write(msg + "\n"); logfile.flush()
    
    return chain, log_prob, sampler.iteration

# ============================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--init", choices=['map5','prior'], default='map5',
                        help="map5=MAP+/-5%% (rapido); prior=uniforme amplio (lento)")
    parser.add_argument("--n_walkers", type=int, default=64)
    parser.add_argument("--n_steps", type=int, default=20000)
    parser.add_argument("--models", default='0,1,2,3')
    args = parser.parse_args()
    
    print(f"="*72)
    print(f"TEST 2 V3 MCMC")
    print(f"  init mode: {args.init}")
    print(f"  walkers x steps: {args.n_walkers} x {args.n_steps}")
    print(f"  models: M{args.models}")
    print(f"="*72)
    
    # Cargar datos
    df, sub_inst_names = load_homogeneous_data()
    n_inst = len(sub_inst_names)
    t = df['bjd'].values
    rv = df['rv_resid_c'].values
    rv_err = np.maximum(df['rv_err'].values, 0.1)
    inst_idx = df['inst_idx'].values
    print(f"\n  Data: N={len(t)}, n_inst={n_inst}, baseline={(t.max()-t.min())/365.25:.2f} yr")
    
    # Cargar MAP V3
    with open(MAP_JSON) as f:
        map_data = json.load(f)
    
    # Configuracion modelos
    model_configs = {
        '0': ('M0', False, False),
        '1': ('M1', False, True),
        '2': ('M2', True,  False),
        '3': ('M3', True,  True),
    }
    
    # Logfile
    log_path = f"{OUTDIR}/test2_mcmc_v3_progress.txt"
    logfile = open(log_path, 'a')
    logfile.write(f"\n\n{'='*72}\n")
    logfile.write(f"START: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    logfile.write(f"init={args.init} walkers={args.n_walkers} steps={args.n_steps}\n")
    logfile.flush()
    
    # Run cada modelo
    results = {}
    for m_id in args.models.split(','):
        model_name, use_planet, use_ring = model_configs[m_id]
        p_map = get_map_for_model(map_data, model_name, n_inst) if args.init == 'map5' else None
        
        chain, log_prob, n_iter = run_mcmc_model(
            t, rv, rv_err, inst_idx, n_inst,
            use_planet, use_ring, model_name,
            args.n_walkers, args.n_steps, args.init, p_map,
            logfile=logfile,
        )
        
        # Guardar cadena
        np.save(f"{OUTDIR}/test2_mcmc_v3_{model_name}.npy", chain)
        np.save(f"{OUTDIR}/test2_mcmc_v3_{model_name}_logprob.npy", log_prob)
        
        # Estadisticas marginales
        if use_planet:
            K_med = np.median(chain[:, 0])
            K_lo, K_hi = np.percentile(chain[:, 0], [16, 84])
            P_med = np.median(chain[:, 1])
            P_lo, P_hi = np.percentile(chain[:, 1], [16, 84])
            e_med = np.median(chain[:, 2])
            e_lo, e_hi = np.percentile(chain[:, 2], [16, 84])
            
            M_samples = msini_from_K_earth(chain[:, 0], chain[:, 1], 0.972, chain[:, 2])
            M_med = np.median(M_samples)
            M_lo, M_hi = np.percentile(M_samples, [16, 84])
            
            K_sigma = K_med / ((K_hi - K_lo) / 2)
            
            print(f"\n  [{model_name}] Posteriors Kepleriano:")
            print(f"    K = {K_med:.3f} (+{K_hi-K_med:.3f} / -{K_med-K_lo:.3f}) m/s")
            print(f"    P = {P_med:.3f} (+{P_hi-P_med:.3f} / -{P_med-P_lo:.3f}) d")
            print(f"    e = {e_med:.3f} (+{e_hi-e_med:.3f} / -{e_med-e_lo:.3f})")
            print(f"    M = {M_med:.2f} (+{M_hi-M_med:.2f} / -{M_med-M_lo:.2f}) M_Earth")
            print(f"    K_sigma = {K_sigma:.2f}")
            
            results[model_name] = {
                'K_median': float(K_med), 'K_16': float(K_lo), 'K_84': float(K_hi),
                'P_median': float(P_med), 'P_16': float(P_lo), 'P_84': float(P_hi),
                'e_median': float(e_med), 'e_16': float(e_lo), 'e_84': float(e_hi),
                'M_median': float(M_med), 'M_16': float(M_lo), 'M_84': float(M_hi),
                'K_sigma': float(K_sigma),
                'log_prob_max': float(np.max(log_prob)),
                'log_prob_median': float(np.median(log_prob)),
                'n_samples': int(len(chain)),
                'n_iter': int(n_iter),
            }
        else:
            results[model_name] = {
                'log_prob_max': float(np.max(log_prob)),
                'log_prob_median': float(np.median(log_prob)),
                'n_samples': int(len(chain)),
                'n_iter': int(n_iter),
            }
        
        if use_ring:
            A_r_med = np.median(chain[:, 5])
            A_r_lo, A_r_hi = np.percentile(chain[:, 5], [16, 84])
            P_r_med = np.median(chain[:, 6])
            results[model_name]['A_ring_median'] = float(A_r_med)
            results[model_name]['A_ring_16'] = float(A_r_lo)
            results[model_name]['A_ring_84'] = float(A_r_hi)
            results[model_name]['P_ring_median'] = float(P_r_med)
            print(f"    A_ring = {A_r_med:.3f} (+{A_r_hi-A_r_med:.3f} / -{A_r_med-A_r_lo:.3f}) m/s")
    
    # Bayesian evidence comparativa via maxima log_prob (proxy)
    print(f"\n{'='*72}")
    print(f"COMPARACION DE EVIDENCIA (maxima log_prob como proxy)")
    print(f"{'='*72}")
    if 'M0' in results:
        ref = results['M0']['log_prob_max']
        for name in ['M0','M1','M2','M3']:
            if name in results:
                lp = results[name]['log_prob_max']
                delta = lp - ref
                log10bf = delta / np.log(10)
                print(f"  {name}: log_prob_max={lp:.3f}  Δ vs M0={delta:+.3f}  log10(BF)={log10bf:+.3f}")
    
    # Guardar summary
    summary = {
        'config': {
            'init_mode': args.init,
            'n_walkers': args.n_walkers,
            'n_steps': args.n_steps,
            'models': args.models,
            'priors': {
                'K': '[0, 20] m/s',
                'P': '[217, 247] d',
                'e': '[0, 0.7] uniforme (AMPLIO)',
                'P_ring': '[400, 500] d',
            },
        },
        'data': {
            'N': int(len(t)),
            'n_inst': n_inst,
            'sub_instruments': sub_inst_names,
            'baseline_yr': float((t.max()-t.min())/365.25),
        },
        'results': results,
    }
    with open(f"{OUTDIR}/test2_mcmc_v3_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n[OK] {OUTDIR}/test2_mcmc_v3_summary.json")
    
    logfile.write(f"END: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    logfile.close()

if __name__ == "__main__":
    main()
