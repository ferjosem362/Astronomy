"""
Detection limit map estilo Zhao 2017 Figura 8.

Inject-recovery sobre grilla densa (P x M x phase):
  - 30 periodos log-espaciados entre 1.5 y 2000 d
  - 12 masas log-espaciadas entre 0.1 y 30 M_Earth
  - 20 fases aleatorias por celda
Total: 7200 GLS calculations

Persiste resultados parciales cada N celdas para resilencia.

Output:
  data_interim/qc/detection_limits_grid.csv  (resultados crudos)
  data_interim/qc/detection_limits_map.png   (figura paper-ready)
"""

import os
import time
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

INPUT  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"
GRID_CSV = f"{OUTDIR}/detection_limits_grid.csv"
PLOT_PNG = f"{OUTDIR}/detection_limits_map.png"

M_STAR_SUN = 0.972
HZ_INNER_DAYS = 75   # ~0.5 AU para Alpha Cen B
HZ_OUTER_DAYS = 350  # ~1.0 AU

# Grilla
N_PERIODS = 30
N_MASSES  = 12
N_PHASES  = 20

P_GRID = np.logspace(np.log10(1.5), np.log10(2000), N_PERIODS)
M_GRID = np.logspace(np.log10(0.1), np.log10(30), N_MASSES)

# GLS frecuencias - reducimos resolucion para acelerar
P_min, P_max = 1.0, 5000.0
nf = 25000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods_grid_lookup = 1.0 / freqs

# FAP threshold para recuperacion
FAP_THRESH = 0.01
# Tolerancia de periodo recuperado (% del injectado)
P_TOL_REL = 0.05

# Zhao 2017 Tabla 4 Figura 8 - extracted points
ZHAO_2017 = {
    'P_d':   [3, 5, 10, 30, 53, 100, 200, 500, 1000],
    'M_min': [3.0, 3.0, 3.0, 3.5, 3.5, 3.5, 5.0, 10, 50]
}

# Cargar
df = pd.read_csv(INPUT)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
df = df.sort_values('bjd').reset_index(drop=True)

t = df['bjd'].values
y_orig = df['rv_resid_centered'].values
dy = np.maximum(df['rv_err'].values, 0.1)
N_pts = len(t)

print(f"[INFO] Master clean: N={N_pts}, "
      f"baseline={(t.max()-t.min())/365.25:.2f} yr")
print(f"[INFO] Grid: {N_PERIODS} P x {N_MASSES} M x {N_PHASES} phases = "
      f"{N_PERIODS*N_MASSES*N_PHASES} GLS")
print(f"[INFO] Persistencia parcial: {GRID_CSV}")

# K para una masa dada, P en dias, e=0, sin(i)=1
def K_from_M(M_earth, P_days, M_star=M_STAR_SUN):
    P_yr = P_days / 365.25
    return 0.6395 * M_earth * M_star**(-2/3) * P_yr**(-1/3)

# Pre-calcular GLS baseline (sin injection)
print("\n[INFO] Calculando GLS baseline (sin injection)...")
ls_base = LombScargle(t, y_orig, dy, normalization='standard')
power_base = ls_base.power(freqs)
print(f"[INFO] Baseline OK")

# ============================================================
# Resume from existing CSV if exists
# ============================================================
if os.path.exists(GRID_CSV):
    df_existing = pd.read_csv(GRID_CSV)
    done_set = set(zip(df_existing['P_inj'].round(4),
                       df_existing['M_inj'].round(4),
                       df_existing['phase_idx']))
    print(f"\n[RESUME] Encontrado {GRID_CSV} con {len(df_existing)} celdas ya hechas")
    results = df_existing.to_dict('records')
else:
    results = []
    done_set = set()

# ============================================================
# Loop principal
# ============================================================
np.random.seed(42)
phases = np.random.uniform(0, 2*np.pi, N_PHASES)

t_start = time.time()
last_save = t_start
save_every_sec = 60  # guardar cada minuto

total_cells = N_PERIODS * N_MASSES * N_PHASES
done = len(results)

print(f"\n[INFO] Comenzando loop principal ({total_cells - done} celdas restantes)...")
print(f"{'idx':>5} {'P':>8} {'M':>6} {'K':>7} {'P_phase':>4} "
      f"{'P_rec':>8} {'pow_rec':>8} {'FAP':>10} {'rec':>4}")
print("-" * 75)

for i_p, P_inj in enumerate(P_GRID):
    for i_m, M_inj in enumerate(M_GRID):
        K_inj = K_from_M(M_inj, P_inj)
        
        # Banda de evaluacion: ±5% del periodo inyectado
        p_lo_eval = P_inj * (1 - P_TOL_REL)
        p_hi_eval = P_inj * (1 + P_TOL_REL)
        in_band = (periods_grid_lookup >= p_lo_eval) & (periods_grid_lookup <= p_hi_eval)
        baseline_pow = np.max(power_base[in_band])
        
        for i_ph, phase in enumerate(phases):
            key = (round(P_inj, 4), round(M_inj, 4), i_ph)
            if key in done_set:
                continue
            
            # Inyectar senial
            signal = K_inj * np.sin(2*np.pi*t/P_inj + phase)
            y_inj = y_orig + signal
            
            # GLS
            ls = LombScargle(t, y_inj, dy, normalization='standard')
            power = ls.power(freqs)
            
            # Pico en banda
            i_peak = np.argmax(power * in_band)
            P_rec = periods_grid_lookup[i_peak]
            pow_rec = power[i_peak]
            fap = ls.false_alarm_probability(pow_rec, method='baluev')
            
            # Criterios de recuperacion
            recovered = (fap < FAP_THRESH) and \
                        (abs(P_rec - P_inj) / P_inj < P_TOL_REL) and \
                        (pow_rec > baseline_pow + 0.01)
            
            results.append({
                'P_inj': P_inj,
                'M_inj': M_inj,
                'K_inj': K_inj,
                'phase_idx': i_ph,
                'phase': phase,
                'P_rec': P_rec,
                'power_rec': pow_rec,
                'baseline_pow': baseline_pow,
                'FAP': fap,
                'recovered': recovered,
            })
            done += 1
            
            # Print cada 100 celdas
            if done % 100 == 0:
                rec_str = "Y" if recovered else "n"
                elapsed = time.time() - t_start
                rate = done / elapsed
                eta_min = (total_cells - done) / rate / 60 if rate > 0 else 0
                print(f"{done:>5}/{total_cells} "
                      f"P={P_inj:>6.1f} M={M_inj:>5.2f} K={K_inj:>5.3f} "
                      f"#{i_ph:<2} P_rec={P_rec:>6.2f} pow={pow_rec:.3f} "
                      f"FAP={fap:.1e} {rec_str}  "
                      f"[ETA {eta_min:.0f} min]")
            
            # Guardar cada 60 sec
            if time.time() - last_save > save_every_sec:
                pd.DataFrame(results).to_csv(GRID_CSV, index=False)
                last_save = time.time()

# Save final
pd.DataFrame(results).to_csv(GRID_CSV, index=False)
elapsed = time.time() - t_start
print(f"\n[OK] {GRID_CSV} ({len(results)} celdas, {elapsed/60:.1f} min)")

# ============================================================
# Plot detection limit map
# ============================================================
print("\n[INFO] Generando figura...")

df_res = pd.DataFrame(results)

# Promediar % recuperado por (P, M)
df_agg = df_res.groupby(['P_inj', 'M_inj']).agg(
    pct_rec=('recovered', 'mean')
).reset_index()

# Pivotar a matriz P x M
P_unique = sorted(df_agg['P_inj'].unique())
M_unique = sorted(df_agg['M_inj'].unique())

Z = np.zeros((len(M_unique), len(P_unique)))
for _, r in df_agg.iterrows():
    i_m = M_unique.index(r['M_inj'])
    i_p = P_unique.index(r['P_inj'])
    Z[i_m, i_p] = r['pct_rec'] * 100  # en %

# Encontrar M_min al 99% por periodo
M_min_99 = []
for i_p, P_v in enumerate(P_unique):
    col = Z[:, i_p]
    # Encontrar la masa minima donde >=99%
    idx = np.where(col >= 99)[0]
    if len(idx) > 0:
        # Interpolar entre la celda anterior y esta
        i_first = idx[0]
        if i_first > 0:
            # Interpolacion log
            f1 = col[i_first-1]
            f2 = col[i_first]
            m1 = np.log10(M_unique[i_first-1])
            m2 = np.log10(M_unique[i_first])
            m_99 = m1 + (99 - f1)/(f2 - f1) * (m2 - m1)
            M_min_99.append(10**m_99)
        else:
            M_min_99.append(M_unique[i_first])
    else:
        M_min_99.append(np.nan)

M_min_99 = np.array(M_min_99)
P_arr = np.array(P_unique)

# Encontrar M_min al 90% y 50% tambien
def m_min_at_pct(Z_col, M_unique_log, target_pct):
    idx = np.where(Z_col >= target_pct)[0]
    if len(idx) == 0: return np.nan
    i_first = idx[0]
    if i_first == 0: return 10**M_unique_log[i_first]
    f1, f2 = Z_col[i_first-1], Z_col[i_first]
    m1, m2 = M_unique_log[i_first-1], M_unique_log[i_first]
    if f2 == f1: return 10**m1
    m_t = m1 + (target_pct - f1)/(f2 - f1) * (m2 - m1)
    return 10**m_t

M_unique_log = np.log10(M_unique)
M_min_50 = np.array([m_min_at_pct(Z[:, i], M_unique_log, 50) for i in range(len(P_unique))])
M_min_90 = np.array([m_min_at_pct(Z[:, i], M_unique_log, 90) for i in range(len(P_unique))])

# Figura
fig, axes = plt.subplots(1, 2, figsize=(16, 6.5),
                          gridspec_kw={'width_ratios': [3, 1.2]})

# === Panel izquierdo: heatmap ===
ax = axes[0]
P_mesh, M_mesh = np.meshgrid(P_unique, M_unique)
pcm = ax.pcolormesh(P_mesh, M_mesh, Z, cmap='viridis',
                     shading='auto', vmin=0, vmax=100)

# Contornos
levels = [50, 90, 99]
colors_c = ['white', 'yellow', 'red']
cs = ax.contour(P_mesh, M_mesh, Z, levels=levels,
                colors=colors_c, linewidths=2)
ax.clabel(cs, inline=True, fontsize=9, fmt='%d%%')

# Linea M_min al 99% (este trabajo)
mask_valid = ~np.isnan(M_min_99)
ax.plot(P_arr[mask_valid], M_min_99[mask_valid], 'r-', lw=2.5,
        label='Este trabajo (99% recovery)', zorder=5)

# Overlay Zhao 2017
ax.plot(ZHAO_2017['P_d'], ZHAO_2017['M_min'], 'r--', lw=2,
        label='Zhao+2017 (99% recovery, HARPS+CHIRON)', zorder=4)

# Banda habitable
ax.axvspan(HZ_INNER_DAYS, HZ_OUTER_DAYS, color='green', alpha=0.15,
           label='Habitable zone (~0.5-1.0 AU)')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Period (days)', fontsize=12)
ax.set_ylabel('M·sin(i) (M$_\\oplus$)', fontsize=12)
ax.set_title(f'Detection limits in α Cen B\n'
             f'19.6 yr baseline, N={N_pts} nightly RVs, '
             f'{N_PHASES} injections per (P,M) cell',
             fontsize=11)
ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
ax.set_xlim(1.5, 2000)
ax.set_ylim(0.1, 30)
ax.grid(alpha=0.3, which='both')

cbar = plt.colorbar(pcm, ax=ax, label='Recovery rate (%)')

# === Panel derecho: cortes verticales en P keys ===
ax2 = axes[1]
P_KEYS = [3, 10, 53, 74, 270, 700]
colors_p = plt.cm.plasma(np.linspace(0, 0.85, len(P_KEYS)))

for P_k, c in zip(P_KEYS, colors_p):
    # Encontrar la columna mas cercana
    i_p = np.argmin(np.abs(np.array(P_unique) - P_k))
    P_actual = P_unique[i_p]
    col = Z[:, i_p]
    ax2.plot(M_unique, col, '-o', color=c, ms=4,
             label=f'P={P_actual:.0f} d')

ax2.axhline(99, color='red', ls=':', lw=1, alpha=0.5)
ax2.axhline(50, color='gray', ls=':', lw=1, alpha=0.5)
ax2.set_xscale('log')
ax2.set_xlabel('M·sin(i) (M$_\\oplus$)')
ax2.set_ylabel('Recovery rate (%)')
ax2.set_title('Slices por periodo')
ax2.legend(fontsize=8, loc='lower right')
ax2.set_ylim(-5, 105)
ax2.grid(alpha=0.3, which='both')

fig.tight_layout()
fig.savefig(PLOT_PNG, dpi=140, bbox_inches='tight')
print(f"[OK] {PLOT_PNG}")

# ============================================================
# Tabla M_min al 99% (texto)
# ============================================================
print("\n" + "="*70)
print("TABLA: M·sini detection limits al 99% recovery")
print("="*70)
print(f"{'P (d)':>8}  {'M_min este (M_E)':>20}  {'Zhao 2017 (M_E)':>18}  {'mejora':>10}")
print("-"*70)

for P_z, M_z in zip(ZHAO_2017['P_d'], ZHAO_2017['M_min']):
    # Encontrar P mas cercano en nuestra grid
    i_p = np.argmin(np.abs(np.array(P_unique) - P_z))
    P_actual = P_unique[i_p]
    M_us = M_min_99[i_p]
    if not np.isnan(M_us):
        improvement = M_z / M_us
        print(f"{P_z:>8}  {M_us:>15.2f} (P={P_actual:.1f})  "
              f"{M_z:>18.1f}  {improvement:>8.1f}x")
    else:
        print(f"{P_z:>8}  {'NaN (no detectable)':>20}  {M_z:>18.1f}  {'-':>10}")
print("="*70)
