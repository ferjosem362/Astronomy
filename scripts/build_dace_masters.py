"""
Capa 2 final - Construye masters DACE limpios para Alpha Cen B.

Para cada sub-instrumento:
1. Aplica filtros: is_latest_drs + drs_qc + outliers de RV
2. Genera version 'full' (exposicion individual)
3. Genera version 'nightly' (binned por noche con error propagado)

Para el master global:
- dace_master_full.csv  (todas las epocas filtradas)
- dace_master_full_post2008.csv  (sin pre-2008)
- dace_master_nightly.csv
- dace_master_nightly_post2008.csv

Salida adicional: plot diagnostico antes/despues filtros.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

QCDIR = "data_interim/qc"

# -----------------------------------------------------------------------
# Configuracion de filtros por sub-instrumento
# -----------------------------------------------------------------------
SUB_INSTS = {
    'HARPS03':         {'csv': f'{QCDIR}/dace_HARPS03_HARPS.csv',
                        'rv_window_ms': 2000.0},
    'HARPS15':         {'csv': f'{QCDIR}/dace_HARPS15_HARPS.csv',
                        'rv_window_ms': 1000.0},
    'ESPRESSO18_HR11': {'csv': f'{QCDIR}/dace_ESPRESSO18_SINGLEHR11.csv',
                        'rv_window_ms': 500.0},
    'ESPRESSO18_HR21': {'csv': f'{QCDIR}/dace_ESPRESSO18_SINGLEHR21.csv',
                        'rv_window_ms': 100.0},
    'ESPRESSO18_UHR':  {'csv': f'{QCDIR}/dace_ESPRESSO18_SINGLEUHR.csv',
                        'rv_window_ms': 300.0},
    'ESPRESSO19_HR11': {'csv': f'{QCDIR}/dace_ESPRESSO19_SINGLEHR11.csv',
                        'rv_window_ms': 500.0},
}

# Ventana RV: |rv - median(rv)| <= window
# 2000 m/s para HARPS03 cubre 11 anos de movimiento binario con holgura
# 100 m/s para HR21 porque solo cubre 1 dia (esencialmente RV constante)

# -----------------------------------------------------------------------
# Columnas a propagar al master limpio
# -----------------------------------------------------------------------
KEEP = [
    'bjd', 'rv', 'rv_err',
    'ccf_bispan', 'ccf_bispan_err',
    'ccf_fwhm', 'ccf_fwhm_err',
    'ccf_contrast', 'ccf_contrast_err',
    'spectro_rhk', 'spectro_rhk_err',
    'spectro_halpha', 'spectro_halpha_err',
    'spectro_na', 'spectro_na_err',
    'spectro_sn50',
    'cal_berv', 'cal_drift_rv',
    'texp', 'date_night',
    'sub_inst', 'inst_split', 'dace_inst', 'dace_drs', 'dace_mode',
]

# -----------------------------------------------------------------------
def load_and_filter(name, cfg):
    """Carga CSV, aplica filtros de calidad, devuelve DF limpio."""
    df = pd.read_csv(cfg['csv'])
    n0 = len(df)
    
    # Filtro 1: latest DRS + QC
    df = df[(df['is_latest_drs'] == True) & (df['drs_qc'] == True)].copy()
    n1 = len(df)
    
    # Filtro 2: ventana RV alrededor de mediana (anti-contaminacion A)
    rv_med = df['rv'].median()
    in_window = (df['rv'] - rv_med).abs() <= cfg['rv_window_ms']
    df = df[in_window].copy()
    n2 = len(df)
    
    # Filtro 3: outliers extremos en eRV (instrumento mal calibrado esa noche)
    erv_med = df['rv_err'].median()
    erv_ok = df['rv_err'] < 10 * erv_med
    df = df[erv_ok].copy()
    n3 = len(df)
    
    print(f"\n[{name}]")
    print(f"  Inicial:                 {n0}")
    print(f"  Tras latest_drs + qc:    {n1}  (-{n0-n1})")
    print(f"  Tras ventana RV +/-{cfg['rv_window_ms']:.0f}: {n2}  (-{n1-n2})")
    print(f"  Tras filtro eRV:         {n3}  (-{n2-n3})")
    print(f"  RV mediana: {df['rv'].median():.3f}  std: {df['rv'].std():.3f}")
    print(f"  eRV mediana: {df['rv_err'].median():.4f}")
    
    # Quedarnos con columnas relevantes
    keep_cols = [c for c in KEEP if c in df.columns]
    df = df[keep_cols].copy()
    
    return df


def night_bin(df):
    """
    Binning por noche.
    BJD se promedia (no exacto pero suficiente para nuestros fines).
    RV se promedia con peso 1/eRV^2; eRV se propaga como 1/sqrt(sum(1/eRV^2)).
    Otros indicadores se promedian con peso 1/eRV^2 (mismo peso que RV).
    """
    df = df.copy()
    # Asignar noche (date_night ya viene de DACE; si no, redondear BJD-0.5)
    if 'date_night' in df.columns and df['date_night'].notna().any():
        df['_night'] = df['date_night']
    else:
        df['_night'] = np.floor(df['bjd'] - 0.5).astype(int)
    
    df['_w'] = 1.0 / df['rv_err']**2
    
    grouped = df.groupby(['_night', 'sub_inst'])
    
    rows = []
    for (night, sub), g in grouped:
        w = g['_w'].values
        sw = w.sum()
        n_exp = len(g)
        row = {
            'date_night': night,
            'sub_inst': sub,
            'n_exposures': n_exp,
            'bjd': (g['bjd'] * w).sum() / sw,
            'rv': (g['rv'] * w).sum() / sw,
            'rv_err': 1.0 / np.sqrt(sw),
        }
        # propagar otros indicadores con mismo esquema
        for col in ['ccf_bispan', 'ccf_fwhm', 'ccf_contrast',
                    'spectro_rhk', 'spectro_halpha', 'spectro_na',
                    'spectro_sn50', 'cal_berv', 'cal_drift_rv', 'texp']:
            if col in g.columns and g[col].notna().any():
                row[col] = (g[col].fillna(0) * w).sum() / sw
        # propagar errores en cuadratura (aprox: 1/sqrt(sum(1/err^2)))
        for ecol in ['ccf_bispan_err', 'ccf_fwhm_err', 'ccf_contrast_err',
                     'spectro_rhk_err', 'spectro_halpha_err', 'spectro_na_err']:
            if ecol in g.columns and g[ecol].notna().any():
                w_ind = 1.0 / g[ecol].fillna(g[ecol].median())**2
                row[ecol] = 1.0 / np.sqrt(w_ind.sum())
        # heredar metadatos del primer punto
        for col in ['inst_split', 'dace_inst', 'dace_drs', 'dace_mode']:
            if col in g.columns:
                row[col] = g[col].iloc[0]
        rows.append(row)
    
    return pd.DataFrame(rows).sort_values('bjd').reset_index(drop=True)


# =======================================================================
# MAIN
# =======================================================================

print("="*78)
print("CAPA 2 FINAL: construyendo masters DACE limpios")
print("="*78)

all_full = []
all_nightly = []

for name, cfg in SUB_INSTS.items():
    df_full = load_and_filter(name, cfg)
    if len(df_full) == 0:
        print(f"  [SKIP] {name} sin datos tras filtros")
        continue
    
    out_full = f"{QCDIR}/dace_{name}_full.csv"
    df_full.to_csv(out_full, index=False)
    
    df_night = night_bin(df_full)
    out_night = f"{QCDIR}/dace_{name}_nightly.csv"
    df_night.to_csv(out_night, index=False)
    
    print(f"  -> {out_full}  (N={len(df_full)})")
    print(f"  -> {out_night}  (N={len(df_night)} noches)")
    
    all_full.append(df_full)
    all_nightly.append(df_night)

# -----------------------------------------------------------------------
# Masters globales
# -----------------------------------------------------------------------
print("\n" + "="*78)
print("MASTERS GLOBALES")
print("="*78)

if all_full:
    df_master_full = pd.concat(all_full, ignore_index=True).sort_values('bjd').reset_index(drop=True)
    df_master_full.to_csv(f"{QCDIR}/dace_master_full.csv", index=False)
    print(f"  dace_master_full.csv:    N={len(df_master_full)}")
    
    bjd_2008 = 2454466.5  # 1 ene 2008
    df_post2008 = df_master_full[df_master_full['bjd'] >= bjd_2008].copy()
    df_post2008.to_csv(f"{QCDIR}/dace_master_full_post2008.csv", index=False)
    print(f"  dace_master_full_post2008.csv:    N={len(df_post2008)}")

if all_nightly:
    df_master_night = pd.concat(all_nightly, ignore_index=True).sort_values('bjd').reset_index(drop=True)
    df_master_night.to_csv(f"{QCDIR}/dace_master_nightly.csv", index=False)
    print(f"  dace_master_nightly.csv: N={len(df_master_night)} noches")
    
    bjd_2008 = 2454466.5
    df_night_post = df_master_night[df_master_night['bjd'] >= bjd_2008].copy()
    df_night_post.to_csv(f"{QCDIR}/dace_master_nightly_post2008.csv", index=False)
    print(f"  dace_master_nightly_post2008.csv: N={len(df_night_post)} noches")

# -----------------------------------------------------------------------
# Plot diagnostico
# -----------------------------------------------------------------------
print("\n" + "="*78)
print("PLOT DIAGNOSTICO")
print("="*78)

fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
colors = {'HARPS03': 'tab:blue', 'HARPS15': 'tab:cyan',
          'ESPRESSO18_HR11': 'tab:red', 'ESPRESSO18_HR21': 'tab:orange',
          'ESPRESSO18_UHR': 'tab:pink', 'ESPRESSO19_HR11': 'tab:purple'}

# Panel 1: full
for df in all_full:
    sub = df['sub_inst'].iloc[0]
    axes[0].errorbar(df['bjd'], df['rv'], yerr=df['rv_err'],
                     fmt='.', ms=2, alpha=0.4,
                     color=colors.get(sub, 'gray'),
                     label=f"{sub} (N={len(df)})")
axes[0].set_ylabel('RV (m/s)')
axes[0].set_title('DACE master FULL (todas exposiciones, post-filtros)')
axes[0].legend(fontsize=8, loc='upper left', ncol=3)
axes[0].grid(alpha=0.3)

# Panel 2: nightly
for df in all_nightly:
    sub = df['sub_inst'].iloc[0]
    axes[1].errorbar(df['bjd'], df['rv'], yerr=df['rv_err'],
                     fmt='o', ms=4, alpha=0.6,
                     color=colors.get(sub, 'gray'),
                     label=f"{sub} (N={len(df)})")
axes[1].set_ylabel('RV (m/s)')
axes[1].set_xlabel('BJD')
axes[1].set_title('DACE master NIGHTLY (binned por noche)')
axes[1].legend(fontsize=8, loc='upper left', ncol=3)
axes[1].grid(alpha=0.3)

# Anotar fiber upgrade y cryo
for ax in axes:
    ax.axvline(2457167.5, color='red', ls='--', lw=1, alpha=0.5)
    ax.axvline(2458645.0, color='orange', ls='--', lw=1, alpha=0.5)

fig.tight_layout()
diag_png = f"{QCDIR}/dace_master_diagnostic.png"
fig.savefig(diag_png, dpi=120)
print(f"  -> {diag_png}")
print("\n[OK] Capa 2 cerrada.")
