"""
Inject-recovery dirigido por hipotesis externa:
  - Hipotesis A: planeta a P=53.5 d, K=1.3 m/s (M·sini ~ 1.08 M_Earth)
  - Hipotesis B: planeta a P=262 d, K=0.85 m/s (M·sini ~ 1.17 M_Earth)

Tres tareas:
  1. Buscar picos REALES en bandas estrechas alrededor de 53.5 d y 262 d
  2. Inyectar grilla densa de K en P=53.5 y P=262
  3. Para el caso especifico de hipotesis (K=1.3 a 53.5; K=0.85 a 262):
     - inyectar la senial
     - verificar si el GLS la recupera con FAP < 1%
     - reportar si los datos *originales* contienen ya esa senial (compatibilidad)

Output:
  data_interim/qc/inject_recovery_external.csv
  data_interim/qc/inject_recovery_external_diagnostic.png
  data_interim/qc/zoom_53d_y_262d.png
"""

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"

M_STAR_SUN   = 0.972
M_MARS_EARTH = 0.107

# Bandas de hipotesis
HYPOTHESES = [
    {'P': 53.5, 'K_hyp': 1.3,  'M_hyp_earth': 1.08, 'label': 'H_A_53.5d'},
    {'P': 262.0, 'K_hyp': 0.85, 'M_hyp_earth': 1.17, 'label': 'H_B_262d'},
]

# Grilla densa de K para mapear sensibilidad
K_GRID = [0.1, 0.3, 0.5, 0.7, 0.85, 1.0, 1.3, 1.6, 2.0, 3.0]

# Frecuencias GLS
P_min, P_max = 1.5, 5000.0
nf = 80000  # mas resolucion para zoom
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods_grid = 1.0 / freqs

# Cargar
df = pd.read_csv(INPUT)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
df = df.sort_values('bjd').reset_index(drop=True)
print(f"[INFO] N={len(df)}, baseline={(df['bjd'].max()-df['bjd'].min())/365.25:.2f} yr")

t = df['bjd'].values
y_orig = df['rv_resid_centered'].values
dy = np.maximum(df['rv_err'].values, 0.1)

# ============================================================
# FASE 1: buscar picos REALES en bandas estrechas
# ============================================================
print("\n" + "="*70)
print("FASE 1: BUSQUEDA DIRIGIDA EN BANDAS DE HIPOTESIS")
print("="*70)

ls_orig = LombScargle(t, y_orig, dy, normalization='standard')
power_orig = ls_orig.power(freqs)

for h in HYPOTHESES:
    P_hyp = h['P']
    p_lo = P_hyp - 5.0  # banda ancha para inspeccion
    p_hi = P_hyp + 5.0
    in_band = (periods_grid >= p_lo) & (periods_grid <= p_hi)
    n_pts = in_band.sum()
    
    if n_pts == 0:
        print(f"\n[{h['label']}] sin puntos en banda")
        continue
    
    pwr_band = power_orig[in_band]
    p_band = periods_grid[in_band]
    i_max = np.argmax(pwr_band)
    P_max_band = p_band[i_max]
    pwr_max_band = pwr_band[i_max]
    fap_max = ls_orig.false_alarm_probability(pwr_max_band, method='baluev')
    
    # Comparar con power esperado de la senial inyectada
    # power_expected ~ (K^2 N) / (4 sigma^2)  para senial perfectamente recuperada
    sigma_med = np.median(dy)
    pwr_expected = (h['K_hyp']**2 * len(t)) / (4 * sigma_med**2)
    # Pero LS standard normaliza diferente, asi que esto es solo orientativo
    
    print(f"\n[{h['label']}] P_hyp={P_hyp} d, K_hyp={h['K_hyp']} m/s, M={h['M_hyp_earth']} M_E")
    print(f"  Banda inspeccionada: {p_lo:.1f} – {p_hi:.1f} d")
    print(f"  Pico mas alto en banda (DATOS REALES): P={P_max_band:.3f} d, "
          f"power={pwr_max_band:.4f}, FAP={fap_max:.2e}")
    if fap_max < 0.01:
        print(f"  ** Hay senial significativa real en esta banda **")
    elif fap_max < 0.10:
        print(f"  * Hay senial marginal en esta banda *")
    else:
        print(f"  - No hay senial significativa real en esta banda")
    
    h['P_real_peak'] = P_max_band
    h['power_real_peak'] = pwr_max_band
    h['fap_real_peak'] = fap_max

# ============================================================
# FASE 2: inject-recovery con grilla densa
# ============================================================
print("\n" + "="*70)
print("FASE 2: INJECT-RECOVERY GRID")
print("="*70)

results = []

for h in HYPOTHESES:
    P_inj = h['P']
    label = h['label']
    
    # Banda estrecha (+/-3% del periodo) para evaluar recuperacion
    p_lo_eval = P_inj * 0.97
    p_hi_eval = P_inj * 1.03
    in_band_eval = (periods_grid >= p_lo_eval) & (periods_grid <= p_hi_eval)
    
    print(f"\n[{label}] P_inj={P_inj} d")
    print(f"  Power baseline (sin inj) en banda: "
          f"{np.max(power_orig[in_band_eval]):.4f}")
    print(f"\n  {'K_inj':>8} {'P_peak':>9} {'power':>9} {'excess':>9} "
          f"{'FAP':>10} {'recov':>6}")
    
    baseline_power = np.max(power_orig[in_band_eval])
    
    for K_inj in K_GRID:
        signal = K_inj * np.sin(2*np.pi * t / P_inj)
        y_inj = y_orig + signal
        ls_inj = LombScargle(t, y_inj, dy, normalization='standard')
        power_inj = ls_inj.power(freqs)
        
        i_peak = np.argmax(power_inj * in_band_eval)
        P_peak = periods_grid[i_peak]
        pow_peak = power_inj[i_peak]
        excess = pow_peak - baseline_power
        fap = ls_inj.false_alarm_probability(pow_peak, method='baluev')
        recovered = (excess > 0.02) and (fap < 0.01) and \
                    (abs(P_peak - P_inj) / P_inj < 0.03)
        rec_str = "YES" if recovered else "no"
        
        # Calcular masa
        P_yr = P_inj / 365.25
        M_p_earth = K_inj / 0.6395 * (M_STAR_SUN)**(2/3) * (P_yr)**(1/3)
        
        print(f"  {K_inj:>8.3f} {P_peak:>9.3f} {pow_peak:>9.4f} "
              f"{excess:>+9.4f} {fap:>10.2e} {rec_str:>6}")
        
        results.append({
            'hypothesis': label,
            'P_inj': P_inj,
            'K_inj': K_inj,
            'M_p_earth': M_p_earth,
            'P_peak': P_peak,
            'power': pow_peak,
            'excess': excess,
            'FAP': fap,
            'recovered': recovered,
            'is_hypothesis_K': abs(K_inj - h['K_hyp']) < 0.01,
        })

df_res = pd.DataFrame(results)
df_res.to_csv(f"{OUTDIR}/inject_recovery_external.csv", index=False)
print(f"\n[OK] {OUTDIR}/inject_recovery_external.csv")

# ============================================================
# FASE 3: foco en el caso especifico de la hipotesis
# ============================================================
print("\n" + "="*70)
print("FASE 3: TEST DEL CASO HIPOTETICO ESPECIFICO")
print("="*70)

for h in HYPOTHESES:
    print(f"\n[{h['label']}]")
    print(f"  Hipotesis: P={h['P']} d, K={h['K_hyp']} m/s, M={h['M_hyp_earth']} M_E")
    
    # Datos REALES (sin inyectar): ¿se ve la senial?
    fap_real = h.get('fap_real_peak', np.nan)
    P_real = h.get('P_real_peak', np.nan)
    print(f"  En datos REALES sin inj: P={P_real:.3f} d, "
          f"FAP={fap_real:.2e}")
    if fap_real < 0.01:
        print(f"  ↑ Senial real significativa detectada en esa banda")
    else:
        print(f"  ↑ NO hay senial significativa real en esa banda")
    
    # Caso inyectado de la hipotesis
    sub = df_res[(df_res['hypothesis'] == h['label']) & df_res['is_hypothesis_K']]
    if len(sub) > 0:
        r = sub.iloc[0]
        print(f"  Inyectando K_hyp={r['K_inj']} m/s en datos REALES:")
        print(f"    -> P_peak recuperado={r['P_peak']:.3f} d, FAP={r['FAP']:.2e}, "
              f"recovered={r['recovered']}")
        
        # Conclusion logica
        if r['recovered'] and fap_real < 0.01:
            print(f"  CONCLUSION: La hipotesis ES COMPATIBLE — la senial existe en datos")
        elif r['recovered'] and fap_real >= 0.01:
            print(f"  CONCLUSION: La hipotesis ES DETECTABLE pero NO ESTA en los datos")
            print(f"               -> tension con la hipotesis externa")
        elif not r['recovered']:
            print(f"  CONCLUSION: La hipotesis NO ES DETECTABLE con este pipeline")

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Top: zoom espectral en bandas
for ax, h in zip(axes[0], HYPOTHESES):
    P_hyp = h['P']
    p_lo = P_hyp - 8
    p_hi = P_hyp + 8
    mask = (periods_grid >= p_lo) & (periods_grid <= p_hi)
    
    ax.plot(periods_grid[mask], power_orig[mask], 'k-', lw=1, label='Datos reales')
    
    # FAP levels
    fap_thresh = ls_orig.false_alarm_level([0.10, 0.01, 0.001], method='baluev')
    ax.axhline(fap_thresh[2], color='red', ls='--', lw=0.6, label='FAP 0.1%')
    ax.axhline(fap_thresh[1], color='orange', ls='--', lw=0.6, label='FAP 1%')
    
    ax.axvline(P_hyp, color='green', ls=':', lw=2, label=f'P_hyp={P_hyp}')
    ax.set_xlabel('Period (d)')
    ax.set_ylabel('GLS power')
    ax.set_title(f'Zoom espectral: {h["label"]}')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

# Bottom: recuperacion vs K (como antes pero por banda)
for ax, h in zip(axes[1], HYPOTHESES):
    label = h['label']
    sub = df_res[df_res['hypothesis'] == label]
    
    rec = sub[sub['recovered']]
    not_rec = sub[~sub['recovered']]
    
    ax.scatter(rec['K_inj'], rec['excess'], s=80, color='tab:green',
               label='Recovered', zorder=3, edgecolor='black')
    ax.scatter(not_rec['K_inj'], not_rec['excess'], s=80, color='tab:red',
               label='Not recovered', zorder=3, edgecolor='black', alpha=0.7)
    
    # Marcar K hipotesis
    ax.axvline(h['K_hyp'], color='green', ls=':', lw=2,
               label=f"K_hyp={h['K_hyp']} m/s")
    ax.axhline(0.02, color='gray', ls='--', lw=1, alpha=0.5)
    
    for _, r in sub.iterrows():
        ax.annotate(f"{r['M_p_earth']:.2f}M⊕",
                    xy=(r['K_inj'], r['excess']),
                    xytext=(5, 5), textcoords='offset points', fontsize=7)
    
    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=0.001)
    ax.set_xlabel('K_inj (m/s)')
    ax.set_ylabel('Power excess')
    ax.set_title(f'Recuperabilidad: {label}')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3, which='both')

fig.tight_layout()
out = f"{OUTDIR}/inject_recovery_external_diagnostic.png"
fig.savefig(out, dpi=120)
print(f"\n[OK] {out}")
