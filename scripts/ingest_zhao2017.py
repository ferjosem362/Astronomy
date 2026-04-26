"""
Ingesta de la Tabla 2 de Zhao 2017 al esquema DACE.

Genera dos sub-instrumentos para Alpha Cen B:
  - ES_zhao2017     (228 puntos para A, 133 para B, 2008-2010)
  - CHIRON_zhao2017 (63 puntos para A, 108 para B,  2011-2013)

Salida: data_interim/qc/zhao2017_alphacenb_full.csv
        data_interim/qc/zhao2017_alphacenb_nightly.csv
"""

import os
import sys
import numpy as np
import pandas as pd

INPUT  = "data_raw/zhao2017_table2.txt"
QCDIR  = "data_interim/qc"
TARGET_STAR = "B"   # nos interesa Alpha Cen B

if not os.path.exists(INPUT):
    sys.exit(f"[ERROR] No encuentro {INPUT}")

# Leer tabla
rows = []
with open(INPUT) as f:
    for line in f:
        parts = line.split()
        if len(parts) != 5:
            continue
        star, jd2440000, rv, err, source = parts
        rows.append({
            'star': star,
            'bjd': float(jd2440000) + 2440000.0,
            'rv': float(rv),
            'rv_err': float(err),
            'instrument': source,
        })

df = pd.DataFrame(rows)
print(f"[INFO] Tabla 2 Zhao 2017 leida: {len(df)} filas totales")
print(f"[INFO] Por estrella e instrumento:")
print(df.groupby(['star', 'instrument']).size())

# Filtrar solo Alpha Cen B
df_b = df[df['star'] == TARGET_STAR].copy()
print(f"\n[INFO] Filtrando solo Alpha Cen {TARGET_STAR}: {len(df_b)} epocas")

# Mapear al esquema DACE
df_b['sub_inst']   = df_b['instrument'].map({'ES': 'ES_zhao2017', 'CHIRON': 'CHIRON_zhao2017'})
df_b['inst_split'] = df_b['sub_inst']
df_b['dace_inst']  = df_b['instrument']
df_b['dace_drs']   = 'Zhao2017'
df_b['dace_mode']  = 'binned_nightly'

# Stats por sub-instrumento
print("\n" + "="*60)
print("ESTADISTICAS POR SUB-INSTRUMENTO")
print("="*60)
for sub in df_b['sub_inst'].unique():
    sub_df = df_b[df_b['sub_inst'] == sub]
    print(f"\n[{sub}]")
    print(f"  N = {len(sub_df)}")
    print(f"  BJD: {sub_df['bjd'].min():.1f} -> {sub_df['bjd'].max():.1f}")
    print(f"  Baseline: {sub_df['bjd'].max() - sub_df['bjd'].min():.1f} d "
          f"= {(sub_df['bjd'].max() - sub_df['bjd'].min())/365.25:.2f} anos")
    print(f"  RV mediana: {sub_df['rv'].median():.2f} m/s")
    print(f"  RV std: {sub_df['rv'].std():.2f} m/s")
    print(f"  RV min/max: {sub_df['rv'].min():.2f} / {sub_df['rv'].max():.2f}")
    print(f"  eRV mediana: {sub_df['rv_err'].median():.3f} m/s")

# Guardar (estos datos YA estan binned por noche en Zhao 2017)
out_full = f"{QCDIR}/zhao2017_alphacenb_full.csv"
df_b = df_b[['bjd', 'rv', 'rv_err', 'sub_inst', 'inst_split',
             'dace_inst', 'dace_drs', 'dace_mode']]
df_b.to_csv(out_full, index=False)
print(f"\n[OK] {out_full}: N={len(df_b)}")

# Como Zhao 2017 ya binnea por noche, el "nightly" es identico al "full"
out_night = f"{QCDIR}/zhao2017_alphacenb_nightly.csv"
df_b.to_csv(out_night, index=False)
print(f"[OK] {out_night}: N={len(df_b)} (identico al full porque Zhao ya binned)")

# Diagnostico visual
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 5))
for sub, color in [('ES_zhao2017', 'tab:green'), ('CHIRON_zhao2017', 'tab:olive')]:
    m = df_b['sub_inst'] == sub
    ax.errorbar(df_b.loc[m, 'bjd'], df_b.loc[m, 'rv'],
                yerr=df_b.loc[m, 'rv_err'],
                fmt='o', ms=4, alpha=0.7, color=color,
                label=f"{sub} (N={m.sum()})")
ax.set_xlabel('BJD')
ax.set_ylabel('RV relativa (m/s)')
ax.set_title('Zhao 2017 Tabla 2 - Alpha Cen B')
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
out_png = f"{QCDIR}/zhao2017_alphacenb_diagnostic.png"
fig.savefig(out_png, dpi=120)
print(f"[OK] {out_png}")
