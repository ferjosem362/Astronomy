"""
Replot del diagnostico DACE con colores correctos por sub-instrumento.
Solo lee los CSV ya generados, no recalcula nada.
"""
import os
import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QCDIR = "data_interim/qc"

# Mapping correcto de nombres de sub-instrumento a colores
COLORS = {
    'HARPS03':         'tab:blue',
    'HARPS15':         'tab:cyan',
    'ESPRESSO18_HR11': 'tab:red',
    'ESPRESSO18_HR21': 'tab:orange',
    'ESPRESSO18_UHR':  'tab:pink',
    'ESPRESSO19_HR11': 'tab:purple',
}

# Cargar masters
df_full  = pd.read_csv(f"{QCDIR}/dace_master_full.csv")
df_night = pd.read_csv(f"{QCDIR}/dace_master_nightly.csv")

print(f"Full master:    {len(df_full)} epocas")
print(f"Nightly master: {len(df_night)} noches")
print(f"Sub_inst values: {df_full['sub_inst'].unique()}")

fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)

# Panel 1: full
for sub, color in COLORS.items():
    m = df_full['sub_inst'] == sub
    if m.sum() == 0:
        continue
    axes[0].errorbar(df_full.loc[m, 'bjd'], df_full.loc[m, 'rv'],
                     yerr=df_full.loc[m, 'rv_err'],
                     fmt='.', ms=3, alpha=0.6, color=color,
                     label=f"{sub} (N={m.sum()})")
axes[0].set_ylabel('RV (m/s)')
axes[0].set_title('DACE master FULL (todas exposiciones, post-filtros)')
axes[0].legend(fontsize=9, loc='upper right', ncol=2)
axes[0].grid(alpha=0.3)

# Panel 2: nightly
for sub, color in COLORS.items():
    m = df_night['sub_inst'] == sub
    if m.sum() == 0:
        continue
    axes[1].errorbar(df_night.loc[m, 'bjd'], df_night.loc[m, 'rv'],
                     yerr=df_night.loc[m, 'rv_err'],
                     fmt='o', ms=5, alpha=0.7, color=color,
                     label=f"{sub} (N={m.sum()} noches)")
axes[1].set_ylabel('RV (m/s)')
axes[1].set_xlabel('BJD')
axes[1].set_title('DACE master NIGHTLY (binned por noche)')
axes[1].legend(fontsize=9, loc='upper right', ncol=2)
axes[1].grid(alpha=0.3)

# Marcadores de upgrades
for ax in axes:
    ax.axvline(2457167.5, color='red', ls='--', lw=1, alpha=0.5,
               label='HARPS fiber upgrade')
    ax.axvline(2458645.0, color='orange', ls='--', lw=1, alpha=0.5,
               label='ESPRESSO cryo upgrade')

fig.tight_layout()
out = f"{QCDIR}/dace_master_diagnostic_v2.png"
fig.savefig(out, dpi=120)
print(f"[OK] {out}")
