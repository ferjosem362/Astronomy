"""
Filtro robusto de outliers post-Pourbaix2016.

Para cada sub-instrumento, calcula MAD del residuo centrado y rechaza
puntos con |resid| > 5 * 1.4826 * MAD.

Salida:
  data_interim/qc/master_v2_nightly_residuals_clean.csv
  data_interim/qc/outliers_filter_diagnostic.png
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT  = "data_interim/qc/master_v2_nightly_residuals.csv"
OUTPUT = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
DIAG   = "data_interim/qc/outliers_filter_diagnostic.png"

NSIGMA = 5.0   # umbral robusto

df = pd.read_csv(INPUT)
print(f"[INFO] Cargado {INPUT}: N={len(df)}")

# Centrar por sub-instrumento
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())

# MAD robusto por sub-instrumento
def robust_sigma(x):
    return 1.4826 * (x - x.median()).abs().median()

# Calcular sigma robusto y flag de outlier por grupo
df['_sigma_robust'] = df.groupby('sub_inst')['rv_resid_centered'].transform(robust_sigma)
df['is_outlier'] = df['rv_resid_centered'].abs() > NSIGMA * df['_sigma_robust']

# Reporte por sub-instrumento
print("\n" + "="*78)
print(f"FILTRADO DE OUTLIERS (umbral {NSIGMA}-sigma robusto)")
print("="*78)
print(f"{'sub_inst':<25} {'N':>4} {'Nout':>5} {'sigma_rob':>10} "
      f"{'std_pre':>9} {'std_post':>9}")
print("-"*78)

clean_rows = []
for sub in df['sub_inst'].unique():
    g = df[df['sub_inst'] == sub]
    n = len(g)
    n_out = g['is_outlier'].sum()
    sigma_r = g['_sigma_robust'].iloc[0]
    std_pre = g['rv_resid_centered'].std()
    g_clean = g[~g['is_outlier']]
    std_post = g_clean['rv_resid_centered'].std() if len(g_clean) > 0 else float('nan')
    print(f"{sub:<25} {n:>4} {n_out:>5} {sigma_r:>10.2f} "
          f"{std_pre:>9.2f} {std_post:>9.2f}")
    clean_rows.append(g_clean)

df_clean = pd.concat(clean_rows, ignore_index=True).sort_values('bjd').reset_index(drop=True)
print("-"*78)
print(f"TOTAL: N={len(df)} -> N_clean={len(df_clean)}  "
      f"(rechazados {len(df)-len(df_clean)} = {100*(len(df)-len(df_clean))/len(df):.1f}%)")

# Guardar
df_clean = df_clean.drop(columns=['_sigma_robust', 'is_outlier'])
df_clean.to_csv(OUTPUT, index=False)
print(f"\n[OK] {OUTPUT}")

# Plot diagnostico
COLORS = {
    'HARPS03':         'tab:blue',
    'HARPS15':         'tab:cyan',
    'ESPRESSO18_HR11': 'tab:red',
    'ESPRESSO18_HR21': 'tab:orange',
    'ESPRESSO18_UHR':  'tab:pink',
    'ESPRESSO19_HR11': 'tab:purple',
    'ES_zhao2017':     'tab:green',
    'CHIRON_zhao2017': 'tab:olive',
}
ORDER = ['HARPS03', 'ES_zhao2017', 'CHIRON_zhao2017', 'HARPS15',
         'ESPRESSO18_HR11', 'ESPRESSO18_HR21', 'ESPRESSO18_UHR',
         'ESPRESSO19_HR11']

fig, ax = plt.subplots(figsize=(14, 5))
for sub in ORDER:
    g = df[df['sub_inst'] == sub]
    if len(g) == 0:
        continue
    g_in = g[~g['is_outlier']]
    g_out = g[g['is_outlier']]
    ax.errorbar(g_in['bjd'], g_in['rv_resid_centered'], yerr=g_in['rv_err'],
                fmt='.', ms=3, alpha=0.6, color=COLORS[sub], label=f"{sub} (N={len(g_in)})")
    if len(g_out) > 0:
        ax.scatter(g_out['bjd'], g_out['rv_resid_centered'].clip(-200, 200),
                   marker='x', s=40, color='red', alpha=0.8)

ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('BJD')
ax.set_ylabel('rv_resid centrado (m/s)')
ax.set_title(f'Master v2 nightly - residuos post-Pourbaix2016 + filtro outliers ({NSIGMA}-sigma robusto)\n'
             f'Outliers (cruces rojas) clipped a +/-200 para visualizacion')
ax.set_ylim(-100, 100)
ax.legend(fontsize=8, loc='upper right', ncol=2)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(DIAG, dpi=120)
print(f"[OK] {DIAG}")
