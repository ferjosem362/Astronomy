"""
Recompute detection limits con comparacion vs Zhao 2017 CORRECTA.

Cambios respecto a v1:
1. Valores Zhao 2017 extraidos de Fig 6d (red noise HARPS, escenario justo)
2. Tambien Fig 5d (white noise HARPS, optimista) para contexto
3. Suavizado del nuestro con interpolacion logistica por periodo
4. Caveats metodologicos en la figura (p-value vs recovery rate)
"""
import sys, os
sys.path.insert(0, '.')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from utils.rv_mass import msini_from_K_earth

INPUT  = "data_interim/qc/detection_limits_grid.csv"
OUTDIR = "data_interim/qc"
OUT_TABLE = f"{OUTDIR}/detection_limits_table_v2.csv"
OUT_PLOT  = f"{OUTDIR}/detection_limits_map_v3_correct.png"

M_STAR = 0.972

# ============================================================
# Valores Zhao 2017 - lectura visual del PDF
# ============================================================
# Fig 5d: HARPS white noise (escenario optimista, sin actividad)
ZHAO_HARPS_WHITE = {
    'P_d':    [3,   10,  30,  53,  100, 200, 500, 1000],
    'M_min':  [3.0, 4.0, 5.0, 6.0, 6.5, 8.0, 15,  28],
}

# Fig 6d: HARPS red noise (escenario realista, INCLUYE actividad)
# Esta es la comparacion JUSTA porque nuestro analisis tambien incluye actividad
ZHAO_HARPS_RED = {
    'P_d':    [3,   10,  30,  53,  100, 200, 500, 1000],
    'M_min':  [5.0, 7.0, 8.0, 9.0, 11,  15,  25,  40],
}

# Fig 5c: CHIRON white noise (peor que HARPS)
ZHAO_CHIRON_WHITE = {
    'P_d':    [3,   10,  30,  100, 300, 1000],
    'M_min':  [4.5, 6.0, 8.0, 12,  20,  35],
}

# ============================================================
def logistic(K, K50, slope):
    """Funcion logistica para fittear recovery rate vs K."""
    return 1.0 / (1.0 + np.exp(-slope * (K - K50)))

def find_K90_logistic(K_vals, recovery_vals):
    """
    Ajusta logistica a (K, recovery) y devuelve K donde recovery=0.9.
    Mas robusto que el primer-bin-superior-a-0.9.
    """
    try:
        # Inicializar K50 en mediana, slope arbitrario
        p0 = [np.median(K_vals), 2.0]
        bounds = ([0, 0.1], [10, 50])
        popt, _ = curve_fit(logistic, K_vals, recovery_vals,
                            p0=p0, bounds=bounds, maxfev=2000)
        K50, slope = popt
        # K donde logistica = 0.9: K = K50 + log(9)/slope
        K90 = K50 + np.log(9) / slope
        return K90, K50, slope
    except Exception:
        return np.nan, np.nan, np.nan

# ============================================================
print("="*72)
print("DETECTION LIMITS V2 - COMPARACION CORRECTA vs ZHAO 2017")
print("="*72)

if not os.path.exists(INPUT):
    print(f"\n[ERROR] No existe {INPUT}")
    sys.exit(1)

df = pd.read_csv(INPUT)
df = df.rename(columns={'P_inj': 'period_d'})
df['e_inj'] = 0.0
print(f"\n  Loaded: {INPUT} shape={df.shape}")
print(f"  Periodos unicos: {df['period_d'].nunique()}")

# Recalcular M*sini con formula correcta
df['M_earth_corrected'] = df.apply(
    lambda r: msini_from_K_earth(r['K_inj'], r['period_d'],
                                  M_star=M_STAR, e=r['e_inj']),
    axis=1
)

# ============================================================
# Para cada periodo, fittear logistica y extraer K90
# ============================================================
print(f"\n  Ajustando logistica por periodo (K90 robusto)...")

detection_table = []
fit_quality = {'good': 0, 'noisy': 0, 'failed': 0}

for P_val in sorted(df['period_d'].unique()):
    sub = df[df['period_d'] == P_val].sort_values('K_inj')
    
    # Bins de K
    K_bins = np.arange(0, sub['K_inj'].max()+0.5, 0.25)
    K_centers = []
    recovery_rates = []
    n_per_bin = []
    for i in range(len(K_bins)-1):
        in_bin = (sub['K_inj'] >= K_bins[i]) & (sub['K_inj'] < K_bins[i+1])
        if in_bin.sum() >= 3:
            K_centers.append((K_bins[i] + K_bins[i+1])/2)
            recovery_rates.append(sub.loc[in_bin, 'recovered'].mean())
            n_per_bin.append(in_bin.sum())
    
    if len(K_centers) < 4:
        fit_quality['failed'] += 1
        continue
    
    K_centers = np.array(K_centers)
    recovery_rates = np.array(recovery_rates)
    
    # Fitear logistica
    K90, K50, slope = find_K90_logistic(K_centers, recovery_rates)
    
    if np.isnan(K90):
        fit_quality['failed'] += 1
        # Fallback: primer K con recovery >= 0.9
        for i, (K, R) in enumerate(zip(K_centers, recovery_rates)):
            if R >= 0.9:
                K90 = K
                break
        if np.isnan(K90):
            continue
    else:
        # Calidad del fit
        residuals = recovery_rates - logistic(K_centers, K50, slope)
        if np.std(residuals) < 0.1:
            fit_quality['good'] += 1
        else:
            fit_quality['noisy'] += 1
    
    M90 = msini_from_K_earth(K90, P_val, M_STAR, 0.0)
    
    detection_table.append({
        'P_d': P_val,
        'K90_logistic': K90,
        'K50_logistic': K50 if not np.isnan(K50) else np.nan,
        'slope_logistic': slope if not np.isnan(slope) else np.nan,
        'M90_Earth': M90,
        'n_bins_used': len(K_centers),
    })

df_table = pd.DataFrame(detection_table)
df_table = df_table.sort_values('P_d').reset_index(drop=True)
df_table.to_csv(OUT_TABLE, index=False)
print(f"  [OK] {OUT_TABLE}")
print(f"\n  Calidad de ajustes:")
print(f"    good (std_resid<0.1): {fit_quality['good']}")
print(f"    noisy:                {fit_quality['noisy']}")
print(f"    failed:               {fit_quality['failed']}")

# ============================================================
# Tabla resumen
# ============================================================
print(f"\n  TABLA DETECTION LIMITS (logistic K90):")
print(f"  {'P (d)':>8} {'K90 (m/s)':>11} {'M90 (M_E)':>11}  {'vs Zhao_red'}")
print(f"  {'-'*55}")

# Interpolar Zhao red noise para comparar
zhao_red_interp = np.interp(
    df_table['P_d'],
    ZHAO_HARPS_RED['P_d'],
    ZHAO_HARPS_RED['M_min']
)

for i, r in df_table.iterrows():
    if i % 3 == 0 or i == len(df_table)-1:  # mostrar 1 de cada 3 + ultima
        zhao_val = zhao_red_interp[i]
        ratio = r['M90_Earth'] / zhao_val
        flag = ""
        if ratio < 0.5: flag = "MUCHO MEJOR"
        elif ratio < 0.8: flag = "mejor"
        elif ratio < 1.2: flag = "similar"
        elif ratio < 2.0: flag = "peor"
        else: flag = "MUCHO PEOR"
        print(f"  {r['P_d']:>8.1f} {r['K90_logistic']:>11.3f} {r['M90_Earth']:>11.2f}  "
              f"Zhao={zhao_val:.1f}  ratio={ratio:.2f}  {flag}")

# ============================================================
# Plot: comparacion completa
# ============================================================
print(f"\n  Generando figura comparativa...")

fig, ax = plt.subplots(1, 1, figsize=(11, 7))

# Nuestros limites
df_clean = df_table.dropna(subset=['M90_Earth'])
ax.plot(df_clean['P_d'], df_clean['M90_Earth'],
        'o-', color='darkred', lw=2, ms=7, alpha=0.85,
        label='This work (master_v3, 90% recovery, formula corrected)')

# Zhao 2017 - HARPS red noise (comparacion JUSTA)
ax.plot(ZHAO_HARPS_RED['P_d'], ZHAO_HARPS_RED['M_min'],
        's--', color='blue', lw=2, ms=8,
        label='Zhao 2017: HARPS red noise (FAP<1%)')

# Zhao 2017 - HARPS white noise (referencia optimista)
ax.plot(ZHAO_HARPS_WHITE['P_d'], ZHAO_HARPS_WHITE['M_min'],
        '^:', color='green', lw=1.5, ms=6, alpha=0.7,
        label='Zhao 2017: HARPS white noise (optimista)')

# Zhao 2017 - CHIRON white noise (peor de Zhao)
ax.plot(ZHAO_CHIRON_WHITE['P_d'], ZHAO_CHIRON_WHITE['M_min'],
        'v:', color='gray', lw=1.5, ms=6, alpha=0.5,
        label='Zhao 2017: CHIRON white noise')

# Zona habitable de alpha Cen B (de Kopparapu 2013)
ax.axvspan(150, 350, color='green', alpha=0.1, label='HZ (Kopparapu+ 2013)')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Period (days)', fontsize=12)
ax.set_ylabel(r'M$\cdot$sin(i) detection limit (M$_\oplus$)', fontsize=12)
ax.set_title(r'Detection limits comparison - $\alpha$ Cen B (M*sini formula corrected)',
             fontsize=12)
ax.legend(loc='upper left', fontsize=9)
ax.grid(alpha=0.3, which='both')
ax.set_xlim(1.5, 1500)
ax.set_ylim(1, 100)

# Caveats
ax.text(0.99, 0.02,
        'Note: This work uses 90% recovery rate threshold;\n'
        'Zhao uses FAP < 1%. Criterios diferentes pero comparables.\n'
        'Zhao values from visual reading of Figs 5d, 5c, 6d.',
        transform=ax.transAxes, fontsize=8, ha='right', va='bottom',
        bbox=dict(facecolor='wheat', alpha=0.7))

fig.tight_layout()
fig.savefig(OUT_PLOT, dpi=120)
print(f"  [OK] {OUT_PLOT}")

# ============================================================
# Resumen final por bandas
# ============================================================
print(f"\n{'='*72}")
print(f"RESUMEN POR BANDAS DE PERIODO (vs Zhao 2017 RED NOISE)")
print(f"{'='*72}")

bands = {
    'Hot (P<10d)':       (1, 10),
    'Warm (10-50d)':     (10, 50),
    'HZ (50-300d)':      (50, 300),
    'Cool (300-1000d)':  (300, 1000),
    'Extra (>1000d)':    (1000, 5000),
}

for bname, (plo, phi) in bands.items():
    in_band = (df_table['P_d'] >= plo) & (df_table['P_d'] < phi)
    if in_band.sum() == 0: continue
    
    M_band = df_table.loc[in_band, 'M90_Earth']
    P_band = df_table.loc[in_band, 'P_d']
    
    # Zhao en este rango
    P_zhao_band = [p for p in ZHAO_HARPS_RED['P_d'] if plo <= p < phi]
    if P_zhao_band:
        M_zhao_band = [ZHAO_HARPS_RED['M_min'][ZHAO_HARPS_RED['P_d'].index(p)] for p in P_zhao_band]
        zhao_med = np.median(M_zhao_band)
    else:
        zhao_med = np.interp((plo+phi)/2, ZHAO_HARPS_RED['P_d'], ZHAO_HARPS_RED['M_min'])
    
    print(f"  {bname:<20} M_min mediana = {M_band.median():>5.1f} M_E   "
          f"Zhao_red = {zhao_med:>5.1f}   ratio = {M_band.median()/zhao_med:.2f}")
