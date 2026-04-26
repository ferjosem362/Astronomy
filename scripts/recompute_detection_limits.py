import sys, os
sys.path.insert(0, '.')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from utils.rv_mass import msini_from_K_earth

INPUT  = "data_interim/qc/detection_limits_grid.csv"
OUTDIR = "data_interim/qc"

OUT_GRID  = f"{OUTDIR}/detection_limits_grid_corrected.csv"
OUT_TABLE = f"{OUTDIR}/detection_limits_table.csv"
OUT_PLOT  = f"{OUTDIR}/detection_limits_map_v3.png"

M_STAR = 0.972

print("="*72)
print("RECALCULAR DETECTION LIMITS - schema P_inj/M_inj/K_inj/recovered")
print("="*72)

if not os.path.exists(INPUT):
    print(f"\n[ERROR] No existe {INPUT}")
    sys.exit(1)

df = pd.read_csv(INPUT)
print(f"\n  Loaded: {INPUT}")
print(f"  Shape: {df.shape}")
print(f"  Cols: {list(df.columns)}")
print(f"\n  Sample (primeras 5 filas):")
print(df.head().to_string())

# Renombrar a convencion del script
df = df.rename(columns={'P_inj': 'period_d'})
df['e_inj'] = 0.0  # asumido circular en el grid original

# Stats originales
print(f"\n  Stats grid:")
print(f"    P range: {df['period_d'].min():.1f} - {df['period_d'].max():.1f} d")
print(f"    K range: {df['K_inj'].min():.3f} - {df['K_inj'].max():.3f} m/s")
print(f"    Recovery overall: {df['recovered'].mean()*100:.1f}%")
print(f"    N periodos unicos: {df['period_d'].nunique()}")

# Recalcular M*sini con formula correcta
df['M_earth_corrected'] = df.apply(
    lambda r: msini_from_K_earth(r['K_inj'], r['period_d'],
                                  M_star=M_STAR, e=r['e_inj']),
    axis=1
)

# Comparar con M_inj antiguo (que tenia el bug)
df = df.rename(columns={'M_inj': 'M_earth_old'})
df['M_factor'] = df['M_earth_corrected'] / df['M_earth_old']
print(f"\n  Factor de correccion (M_corrected / M_old):")
print(f"    median: {df['M_factor'].median():.3f}")
print(f"    min:    {df['M_factor'].min():.3f}")
print(f"    max:    {df['M_factor'].max():.3f}")
print(f"    Esperado: ~7.149 (constante)")

df.to_csv(OUT_GRID, index=False)
print(f"\n  [OK] {OUT_GRID}")

# Detection limits por periodo: K minimo con recuperacion >= 90%
print(f"\n  Calculando detection limits 90%...")

detection_table = []
for P_val in sorted(df['period_d'].unique()):
    sub = df[df['period_d'] == P_val].sort_values('K_inj')
    K_bins = np.arange(0, sub['K_inj'].max()+0.5, 0.25)
    found = False
    for i in range(len(K_bins)-1):
        in_bin = (sub['K_inj'] >= K_bins[i]) & (sub['K_inj'] < K_bins[i+1])
        if in_bin.sum() < 3: continue
        frac = sub.loc[in_bin, 'recovered'].mean()
        if frac >= 0.90:
            K_min = K_bins[i]
            M_min = msini_from_K_earth(K_min, P_val, M_STAR, 0.0)
            detection_table.append({
                'P_d': P_val,
                'K_min_90': K_min,
                'M_min_90_Earth': M_min,
            })
            found = True
            break
    if not found:
        # Si no llega a 90% ni en el K maximo, reportar el K_max
        K_max = sub['K_inj'].max()
        M_max = msini_from_K_earth(K_max, P_val, M_STAR, 0.0)
        detection_table.append({
            'P_d': P_val,
            'K_min_90': np.nan,  # no se alcanzo 90%
            'M_min_90_Earth': np.nan,
        })

df_table = pd.DataFrame(detection_table)
df_table.to_csv(OUT_TABLE, index=False)
print(f"  [OK] {OUT_TABLE}")

# Mostrar tabla
print(f"\n  TABLA DETECTION LIMITS (recuperacion >= 90%):")
print(f"  {'P (d)':>10} {'K_min (m/s)':>15} {'M_min (M_Earth)':>18}")
print(f"  {'-'*45}")
for _, r in df_table.iterrows():
    if pd.isna(r['K_min_90']):
        print(f"  {r['P_d']:>10.1f} {'no llega 90%':>15} {'-':>18}")
    else:
        print(f"  {r['P_d']:>10.1f} {r['K_min_90']:>15.3f} {r['M_min_90_Earth']:>18.2f}")

# Generar figura
print(f"\n  Generando figura...")

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# Panel izquierdo: recovery rate vs K para varios periodos
P_vals = sorted(df['period_d'].unique())
n_show = min(8, len(P_vals))
step = max(1, len(P_vals)//n_show)

for P_val in P_vals[::step]:
    sub = df[df['period_d'] == P_val].sort_values('K_inj')
    K_bins = np.arange(0, sub['K_inj'].max()+0.5, 0.25)
    recovery_rate = []
    for i in range(len(K_bins)-1):
        in_bin = (sub['K_inj'] >= K_bins[i]) & (sub['K_inj'] < K_bins[i+1])
        if in_bin.sum() > 0:
            recovery_rate.append((K_bins[i], sub.loc[in_bin, 'recovered'].mean()))
    if recovery_rate:
        arr = np.array(recovery_rate)
        axes[0].plot(arr[:,0], arr[:,1], '.-', label=f'P={P_val:.0f}d', alpha=0.7, lw=1.2)

axes[0].axhline(0.9, color='red', ls='--', alpha=0.5, lw=1.5, label='90% threshold')
axes[0].set_xlabel('K injected (m/s)')
axes[0].set_ylabel('Recovery fraction')
axes[0].set_title('Recovery rate vs K injected')
axes[0].legend(fontsize=7, loc='lower right', ncol=2)
axes[0].grid(alpha=0.3)
axes[0].set_xlim(0, 5)
axes[0].set_ylim(-0.05, 1.05)

# Panel derecho: detection limits curve M vs P
df_clean = df_table.dropna(subset=['M_min_90_Earth'])
if len(df_clean) > 0:
    axes[1].plot(df_clean['P_d'], df_clean['M_min_90_Earth'],
                 'o-', color='darkred', lw=2, ms=8,
                 label='90% recovery (this work, formula corrected)')
    
    # Comparar con Zhao 2017
    zhao_P = [3, 10, 53, 200, 1000]
    zhao_M = [3.0, 3.0, 3.5, 5.0, 50.0]
    axes[1].plot(zhao_P, zhao_M, 's--', color='gray', lw=1.5, ms=6,
                 label='Zhao et al. 2017 (published)')
    
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].set_xlabel('Period (days)')
    axes[1].set_ylabel('M sin(i) detection limit (M_Earth)')
    axes[1].set_title('Detection limits comparison (corrected formula)')
    axes[1].legend(loc='upper left')
    axes[1].grid(alpha=0.3, which='both')

fig.suptitle('Detection limits alpha Cen B - master_v3 (M*sini corrected, factor 7.15)', y=1.001)
fig.tight_layout()
fig.savefig(OUT_PLOT, dpi=120)
print(f"  [OK] {OUT_PLOT}")

print(f"\n{'='*72}")
print(f"RESUMEN")
print(f"{'='*72}")
print(f"  Grid corrected:  {OUT_GRID}")
print(f"  Tabla:           {OUT_TABLE}")
print(f"  Plot:            {OUT_PLOT}")
print(f"  Factor median:   {df['M_factor'].median():.3f}x")
