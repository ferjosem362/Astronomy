"""
Refinamiento del detection limit map: anyade celdas extra
concentradas en la banda umbral 30%-99% del CSV existente.

Estrategia:
  1. Lee detection_limits_grid.csv ya calculado (7200 celdas)
  2. Para cada periodo, identifica las masas en banda umbral
     (recovery rate entre 30% y 99%)
  3. Anyade 4 masas intermedias entre cada par adyacente en banda umbral
  4. Recalcula 20 fases por celda nueva
  5. Reescribe CSV completo y rehace figura

Output:
  data_interim/qc/detection_limits_grid.csv (extendido)
  data_interim/qc/detection_limits_map.png  (refinado)
  data_interim/qc/detection_limits_map_v2.png (con nueva linea)
"""
import os, time
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_RV  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR    = "data_interim/qc"
GRID_CSV  = f"{OUTDIR}/detection_limits_grid.csv"
PLOT_PNG  = f"{OUTDIR}/detection_limits_map_v2.png"

M_STAR_SUN = 0.972
N_PHASES = 20
P_TOL_REL = 0.05
FAP_THRESH = 0.01

# Cargar datos RV
df = pd.read_csv(INPUT_RV)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
df = df.sort_values('bjd').reset_index(drop=True)
t = df['bjd'].values
y_orig = df['rv_resid_centered'].values
dy = np.maximum(df['rv_err'].values, 0.1)
N_pts = len(t)

# Frecuencias GLS (mismas que en el script original)
P_min, P_max = 1.0, 5000.0
nf = 25000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods_grid_lookup = 1.0 / freqs

print(f"[INFO] Datos RV cargados: N={N_pts}")

# Cargar grid existente
df_grid = pd.read_csv(GRID_CSV)
print(f"[INFO] Grid existente: {len(df_grid)} celdas")

# Agregar por (P, M)
df_agg = df_grid.groupby(['P_inj', 'M_inj']).agg(
    pct_rec=('recovered', 'mean')
).reset_index()
df_agg['pct_rec'] *= 100  # a %

# Identificar (P, M_lo, M_hi) donde anyadir M intermedias
P_values = sorted(df_agg['P_inj'].unique())
new_cells_to_add = []  # lista de (P_inj, M_new)

print(f"\n[INFO] Identificando banda umbral 30%-99% en cada periodo...")
n_thresh_total = 0
for P_v in P_values:
    sub = df_agg[df_agg['P_inj'] == P_v].sort_values('M_inj')
    masses = sub['M_inj'].values
    pcts = sub['pct_rec'].values
    
    # Identificar transiciones cerca del umbral
    in_thresh = (pcts >= 30) & (pcts <= 99)
    n_thresh = in_thresh.sum()
    n_thresh_total += n_thresh
    
    # Para cada par consecutivo donde uno o ambos esta en banda umbral,
    # anyadir 3 masas intermedias log-espaciadas
    for i in range(len(masses) - 1):
        if in_thresh[i] or in_thresh[i+1]:
            m_lo, m_hi = masses[i], masses[i+1]
            log_m = np.linspace(np.log10(m_lo), np.log10(m_hi), 5)[1:-1]  # 3 puntos intermedios
            for log_m_v in log_m:
                new_cells_to_add.append((P_v, 10**log_m_v))

print(f"[INFO] Celdas en banda umbral 30-99%: {n_thresh_total}")
print(f"[INFO] Nuevas celdas a calcular: {len(new_cells_to_add)} x {N_PHASES} fases = "
      f"{len(new_cells_to_add)*N_PHASES} GLS")

# Pre-calcular GLS baseline
print("\n[INFO] Calculando GLS baseline...")
ls_base = LombScargle(t, y_orig, dy, normalization='standard')
power_base = ls_base.power(freqs)

def K_from_M(M_earth, P_days, M_star=M_STAR_SUN):
    P_yr = P_days / 365.25
    return 0.6395 * M_earth * M_star**(-2/3) * P_yr**(-1/3)

# Loop principal sobre nuevas celdas
np.random.seed(123)  # distinto al original para fases diferentes
phases = np.random.uniform(0, 2*np.pi, N_PHASES)

new_results = []
total = len(new_cells_to_add) * N_PHASES
done = 0
t_start = time.time()
last_save = t_start

print(f"\n[INFO] Comenzando loop ({total} GLS)...")

for cell_idx, (P_inj, M_inj) in enumerate(new_cells_to_add):
    K_inj = K_from_M(M_inj, P_inj)
    p_lo_eval = P_inj * (1 - P_TOL_REL)
    p_hi_eval = P_inj * (1 + P_TOL_REL)
    in_band = (periods_grid_lookup >= p_lo_eval) & (periods_grid_lookup <= p_hi_eval)
    baseline_pow = np.max(power_base[in_band])
    
    for i_ph, phase in enumerate(phases):
        signal = K_inj * np.sin(2*np.pi*t/P_inj + phase)
        y_inj = y_orig + signal
        
        ls = LombScargle(t, y_inj, dy, normalization='standard')
        power = ls.power(freqs)
        
        i_peak = np.argmax(power * in_band)
        P_rec = periods_grid_lookup[i_peak]
        pow_rec = power[i_peak]
        fap = ls.false_alarm_probability(pow_rec, method='baluev')
        
        recovered = (fap < FAP_THRESH) and \
                    (abs(P_rec - P_inj) / P_inj < P_TOL_REL) and \
                    (pow_rec > baseline_pow + 0.01)
        
        new_results.append({
            'P_inj': P_inj, 'M_inj': M_inj, 'K_inj': K_inj,
            'phase_idx': i_ph + 100,  # offset para distinguir de las originales
            'phase': phase,
            'P_rec': P_rec, 'power_rec': pow_rec,
            'baseline_pow': baseline_pow,
            'FAP': fap, 'recovered': recovered,
        })
        done += 1
        
        if done % 100 == 0:
            elapsed = time.time() - t_start
            rate = done / elapsed
            eta_min = (total - done) / rate / 60 if rate > 0 else 0
            print(f"  {done}/{total}  P={P_inj:>6.1f} M={M_inj:>5.2f} "
                  f"K={K_inj:.3f}  [ETA {eta_min:.1f} min]")
        
        if time.time() - last_save > 60:
            df_combined = pd.concat([df_grid, pd.DataFrame(new_results)], ignore_index=True)
            df_combined.to_csv(GRID_CSV, index=False)
            last_save = time.time()

# Save final
df_combined = pd.concat([df_grid, pd.DataFrame(new_results)], ignore_index=True)
df_combined.to_csv(GRID_CSV, index=False)
elapsed = time.time() - t_start
print(f"\n[OK] {GRID_CSV} ({len(df_combined)} celdas totales, refinamiento {elapsed/60:.1f} min)")

# ============================================================
# Replot con grid extendido
# ============================================================
print("\n[INFO] Generando figura refinada...")

# Re-aggregate con nuevos puntos
df_agg2 = df_combined.groupby(['P_inj', 'M_inj']).agg(
    pct_rec=('recovered', 'mean'), n_phases=('recovered', 'count')
).reset_index()
df_agg2['pct_rec'] *= 100

# Para heatmap necesitamos grilla regular - interpolamos
P_unique = sorted(df_agg2['P_inj'].unique())
M_unique_full = sorted(df_agg2['M_inj'].unique())

# Construir matriz, NaN donde no hay punto
Z = np.full((len(M_unique_full), len(P_unique)), np.nan)
for _, r in df_agg2.iterrows():
    i_m = M_unique_full.index(r['M_inj'])
    i_p = P_unique.index(r['P_inj'])
    Z[i_m, i_p] = r['pct_rec']

# Calcular linea M_min al 99%, 90%, 50% por interpolacion
def m_min_at_pct(P_v, target_pct):
    sub = df_agg2[df_agg2['P_inj'] == P_v].sort_values('M_inj')
    if len(sub) < 2:
        return np.nan
    masses = sub['M_inj'].values
    pcts = sub['pct_rec'].values
    # Encontrar primer cruce ascendente del target
    above = pcts >= target_pct
    if not above.any():
        return np.nan
    i_first = np.where(above)[0][0]
    if i_first == 0:
        return masses[0]
    p1, p2 = pcts[i_first-1], pcts[i_first]
    m1, m2 = np.log10(masses[i_first-1]), np.log10(masses[i_first])
    if p2 == p1:
        return 10**m1
    log_m_t = m1 + (target_pct - p1)/(p2 - p1) * (m2 - m1)
    return 10**log_m_t

P_arr = np.array(P_unique)
M_min_99 = np.array([m_min_at_pct(P_v, 99) for P_v in P_unique])
M_min_90 = np.array([m_min_at_pct(P_v, 90) for P_v in P_unique])
M_min_50 = np.array([m_min_at_pct(P_v, 50) for P_v in P_unique])

# Zhao 2017
ZHAO_2017 = {
    'P_d':   [3, 5, 10, 30, 53, 100, 200, 500, 1000],
    'M_min': [3.0, 3.0, 3.0, 3.5, 3.5, 3.5, 5.0, 10, 50]
}

# Plot
fig, axes = plt.subplots(1, 2, figsize=(16, 6.5),
                          gridspec_kw={'width_ratios': [3, 1.2]})

ax = axes[0]
P_mesh, M_mesh = np.meshgrid(P_unique, M_unique_full)
pcm = ax.pcolormesh(P_mesh, M_mesh, Z, cmap='viridis',
                     shading='auto', vmin=0, vmax=100)

# Lineas de detection threshold (suaves)
mask99 = ~np.isnan(M_min_99)
mask90 = ~np.isnan(M_min_90)
mask50 = ~np.isnan(M_min_50)
ax.plot(P_arr[mask99], M_min_99[mask99], 'r-', lw=2.5,
        label='99% recovery (este trabajo)', zorder=5)
ax.plot(P_arr[mask90], M_min_90[mask90], '-', color='orange', lw=1.8,
        label='90% recovery', zorder=4)
ax.plot(P_arr[mask50], M_min_50[mask50], '-', color='white', lw=1.5,
        label='50% recovery', zorder=3, alpha=0.7)

# Zhao overlay
ax.plot(ZHAO_2017['P_d'], ZHAO_2017['M_min'], 'r--', lw=2,
        label='Zhao+2017 (99% rec)', zorder=4)

# Habitable zone
ax.axvspan(75, 350, color='green', alpha=0.15, label='Habitable zone')

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('Period (days)', fontsize=12)
ax.set_ylabel('M·sin(i) (M$_\\oplus$)', fontsize=12)
ax.set_title(f'Detection limits in α Cen B (refined grid)\n'
             f'19.6 yr, N={N_pts}, {len(df_combined)} GLS injections',
             fontsize=11)
ax.legend(loc='upper left', fontsize=9, framealpha=0.9)
ax.set_xlim(1.5, 2000); ax.set_ylim(0.1, 30)
ax.grid(alpha=0.3, which='both')
plt.colorbar(pcm, ax=ax, label='Recovery rate (%)')

# Slice panel
ax2 = axes[1]
P_KEYS = [3, 10, 53, 74, 270, 700]
colors_p = plt.cm.plasma(np.linspace(0, 0.85, len(P_KEYS)))

for P_k, c in zip(P_KEYS, colors_p):
    i_p = np.argmin(np.abs(P_arr - P_k))
    P_actual = P_unique[i_p]
    sub = df_agg2[df_agg2['P_inj'] == P_actual].sort_values('M_inj')
    ax2.plot(sub['M_inj'], sub['pct_rec'], '-o', color=c, ms=3,
             label=f'P={P_actual:.0f} d')

ax2.axhline(99, color='red', ls=':', lw=1, alpha=0.5)
ax2.axhline(50, color='gray', ls=':', lw=1, alpha=0.5)
ax2.set_xscale('log')
ax2.set_xlabel('M·sin(i) (M$_\\oplus$)')
ax2.set_ylabel('Recovery rate (%)')
ax2.set_title('Slices por periodo (refinados)')
ax2.legend(fontsize=8, loc='lower right')
ax2.set_ylim(-5, 105)
ax2.grid(alpha=0.3, which='both')

fig.tight_layout()
fig.savefig(PLOT_PNG, dpi=140, bbox_inches='tight')
print(f"[OK] {PLOT_PNG}")

# Tabla actualizada
print("\n" + "="*78)
print("TABLA REFINADA: M·sini detection limits al 99% recovery")
print("="*78)
print(f"{'P (d)':>8}  {'M_min este (M_E)':>22}  {'Zhao 2017 (M_E)':>18}  {'mejora':>10}")
print("-"*78)
for P_z, M_z in zip(ZHAO_2017['P_d'], ZHAO_2017['M_min']):
    i_p = np.argmin(np.abs(P_arr - P_z))
    P_actual = P_unique[i_p]
    M_us = M_min_99[i_p]
    if not np.isnan(M_us):
        improvement = M_z / M_us
        print(f"{P_z:>8}  {M_us:>15.2f} (P={P_actual:.1f})  "
              f"{M_z:>18.1f}  {improvement:>8.1f}x")
    else:
        print(f"{P_z:>8}  {'NaN (insuficiente)':>22}  {M_z:>18.1f}  {'-':>10}")
print("="*78)
