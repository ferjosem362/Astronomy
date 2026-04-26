"""
Master v2 - Cierre de Capa 2.

Combina:
  - DACE (HARPS03, HARPS15, ESPRESSO18 HR11/HR21/UHR, ESPRESSO19 HR11)
  - Zhao 2017 (ES, CHIRON)

Genera 4 masters:
  - master_v2_full.csv             (todas exposiciones DACE + Zhao binned)
  - master_v2_full_post2008.csv    (sin pre-2008 de HARPS03)
  - master_v2_nightly.csv          (DACE binned por noche + Zhao)
  - master_v2_nightly_post2008.csv (sin pre-2008)

Y un diagnostico visual final con colores y leyenda correctos.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QCDIR = "data_interim/qc"
BJD_2008 = 2454466.5

# Mapping para renormalizar nombres de sub_inst a versiones cortas
RENAME = {
    'HARPS03_HARPS':         'HARPS03',
    'HARPS15_HARPS':         'HARPS15',
    'ESPRESSO18_SINGLEHR11': 'ESPRESSO18_HR11',
    'ESPRESSO18_SINGLEHR21': 'ESPRESSO18_HR21',
    'ESPRESSO18_SINGLEUHR':  'ESPRESSO18_UHR',
    'ESPRESSO19_SINGLEHR11': 'ESPRESSO19_HR11',
    # Zhao ya tiene nombres cortos, pero por consistencia
    'ES_zhao2017':           'ES_zhao2017',
    'CHIRON_zhao2017':       'CHIRON_zhao2017',
}

# Colores definitivos (usados en plots y en futuros pasos)
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

# Orden temporal aproximado (para plots y para reportes)
ORDER = ['HARPS03', 'ES_zhao2017', 'CHIRON_zhao2017',
         'HARPS15', 'ESPRESSO18_HR11', 'ESPRESSO18_HR21',
         'ESPRESSO18_UHR', 'ESPRESSO19_HR11']


def load_and_normalize(path):
    df = pd.read_csv(path)
    if 'sub_inst' in df.columns:
        df['sub_inst']   = df['sub_inst'].replace(RENAME)
    if 'inst_split' in df.columns:
        df['inst_split'] = df['inst_split'].replace(RENAME)
    return df


print("="*78)
print("BUILDING MASTER v2 (DACE + Zhao 2017)")
print("="*78)

# -----------------------------------------------------------------
# Cargar DACE
# -----------------------------------------------------------------
df_dace_full  = load_and_normalize(f"{QCDIR}/dace_master_full.csv")
df_dace_night = load_and_normalize(f"{QCDIR}/dace_master_nightly.csv")
print(f"\n[DACE]")
print(f"  full:    {len(df_dace_full)} epocas, sub_inst: {df_dace_full['sub_inst'].unique()}")
print(f"  nightly: {len(df_dace_night)} noches")

# -----------------------------------------------------------------
# Cargar Zhao
# -----------------------------------------------------------------
df_zhao_full = load_and_normalize(f"{QCDIR}/zhao2017_alphacenb_full.csv")
df_zhao_night = df_zhao_full.copy()  # ya viene binned por noche
print(f"\n[Zhao 2017]")
print(f"  N: {len(df_zhao_full)}, sub_inst: {df_zhao_full['sub_inst'].unique()}")

# -----------------------------------------------------------------
# Unificar columnas y concatenar
# -----------------------------------------------------------------
COMMON_COLS = ['bjd', 'rv', 'rv_err', 'sub_inst', 'inst_split']

df_full = pd.concat(
    [df_dace_full[COMMON_COLS + [c for c in df_dace_full.columns if c not in COMMON_COLS]],
     df_zhao_full[COMMON_COLS]],
    ignore_index=True, sort=False
).sort_values('bjd').reset_index(drop=True)

df_night = pd.concat(
    [df_dace_night[COMMON_COLS + [c for c in df_dace_night.columns if c not in COMMON_COLS]],
     df_zhao_night[COMMON_COLS]],
    ignore_index=True, sort=False
).sort_values('bjd').reset_index(drop=True)

# -----------------------------------------------------------------
# Splits temporales
# -----------------------------------------------------------------
df_full_post2008  = df_full[df_full['bjd'] >= BJD_2008].copy()
df_night_post2008 = df_night[df_night['bjd'] >= BJD_2008].copy()

# -----------------------------------------------------------------
# Guardar
# -----------------------------------------------------------------
outputs = {
    'master_v2_full.csv':            df_full,
    'master_v2_full_post2008.csv':   df_full_post2008,
    'master_v2_nightly.csv':         df_night,
    'master_v2_nightly_post2008.csv': df_night_post2008,
}

print("\n" + "="*78)
print("MASTERS GUARDADOS")
print("="*78)
for fname, df in outputs.items():
    out = f"{QCDIR}/{fname}"
    df.to_csv(out, index=False)
    counts = df['sub_inst'].value_counts().to_dict()
    counts_ord = {k: counts.get(k, 0) for k in ORDER if k in counts}
    print(f"\n  {fname}: N={len(df)}")
    for k, v in counts_ord.items():
        print(f"      {k:25s} {v:>6}")

# -----------------------------------------------------------------
# Diagnostico visual final
# -----------------------------------------------------------------
print("\n" + "="*78)
print("PLOT DIAGNOSTICO FINAL")
print("="*78)

fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)

# Panel 1: full
for sub in ORDER:
    m = df_full['sub_inst'] == sub
    if m.sum() == 0:
        continue
    axes[0].errorbar(df_full.loc[m, 'bjd'], df_full.loc[m, 'rv'],
                     yerr=df_full.loc[m, 'rv_err'],
                     fmt='.', ms=3, alpha=0.6, color=COLORS[sub],
                     label=f"{sub} (N={m.sum()})")
axes[0].set_ylabel('RV (m/s)')
axes[0].set_title(f'Master v2 FULL (DACE + Zhao 2017) - N={len(df_full)}')
axes[0].legend(fontsize=8, loc='upper right', ncol=2)
axes[0].grid(alpha=0.3)

# Panel 2: nightly
for sub in ORDER:
    m = df_night['sub_inst'] == sub
    if m.sum() == 0:
        continue
    axes[1].errorbar(df_night.loc[m, 'bjd'], df_night.loc[m, 'rv'],
                     yerr=df_night.loc[m, 'rv_err'],
                     fmt='o', ms=4, alpha=0.7, color=COLORS[sub],
                     label=f"{sub} (N={m.sum()} noches)")
axes[1].set_ylabel('RV (m/s)')
axes[1].set_xlabel('BJD')
axes[1].set_title(f'Master v2 NIGHTLY (DACE binned + Zhao) - N={len(df_night)}')
axes[1].legend(fontsize=8, loc='upper right', ncol=2)
axes[1].grid(alpha=0.3)

# Marcadores
for ax in axes:
    ax.axvline(2457167.5, color='red',    ls='--', lw=1, alpha=0.5)
    ax.axvline(2458645.0, color='orange', ls='--', lw=1, alpha=0.5)
    ax.axvline(BJD_2008,   color='gray',   ls=':',  lw=1, alpha=0.5)

fig.tight_layout()
out_png = f"{QCDIR}/master_v2_diagnostic.png"
fig.savefig(out_png, dpi=120)
print(f"  -> {out_png}")

print("\n" + "="*78)
print("INGESTA CERRADA")
print("="*78)
print(f"\nCobertura temporal:")
print(f"  Total noches (master_v2_nightly):          {len(df_night)}")
print(f"  Tras pre-2008 cut:                         {len(df_night_post2008)}")
print(f"  Baseline total: {df_night['bjd'].max() - df_night['bjd'].min():.0f} dias = "
      f"{(df_night['bjd'].max() - df_night['bjd'].min())/365.25:.1f} anos")
print(f"\nProximo paso sugerido: Fase 3 - periodogramas multi-instrumento.")
