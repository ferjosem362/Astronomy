"""
Integracion de ESPRESSO 2023-2024 expandido (pipeline paralelo Miguel).

Pasos:
  1. Cargar archivo A (rv crudas, 1333 puntos)
  2. Aplicar Pourbaix2016 (mismo modelo que master_v2)
  3. Validar std post-Pourbaix
  4. Filtrar outliers MAD 5sigma
  5. Apilar a nightly (mediana por noche)
  6. Comparar con ESPRESSO19_HR11 actual
  7. Generar diagnostico - NO mergear todavia

Output:
  data_interim/qc/espresso_external_post_pourbaix.csv (per-observation)
  data_interim/qc/espresso_external_nightly.csv       (per-night)
  data_interim/qc/espresso_external_diagnostic.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT  = "data_raw/external_espresso_pipeline_miguel/espresso_2023_2024_curated_hr.csv"
OUTDIR = "data_interim/qc"

# Parametros Pourbaix2016 (mismos que en subtract_binary_pourbaix2016.py)
P_yr_orb = 79.91
e_orb    = 0.524
i_rad    = np.deg2rad(79.32)
omega_B  = np.deg2rad(232.3)
T0_BJD   = 2435349.815
V0_ms    = -22390.0
dVB      = 329.0
K_B      = 5503.97  # m/s

# 1. Cargar
print(f"\n{'='*72}")
print("INTEGRACION ESPRESSO 2023-2024 EXPANDIDO (1333 obs)")
print('='*72)

df = pd.read_csv(INPUT)
print(f"\n[1/6] Cargado {INPUT}")
print(f"      N={len(df)}, columnas: {list(df.columns)}")
print(f"      BJD: {df['time'].min():.0f} - {df['time'].max():.0f} "
      f"({(df['time'].max()-df['time'].min())/365.25:.2f} yr)")
print(f"      Sub-segmentos: {sorted(df['model_inst'].unique())}")
for seg in sorted(df['model_inst'].unique()):
    n = (df['model_inst']==seg).sum()
    print(f"        {seg}: N={n}")

# 2. Aplicar Pourbaix2016
def true_anomaly(t, P_d, e, T0):
    M = 2*np.pi*(t - T0)/P_d % (2*np.pi)
    E = M.copy()
    for _ in range(150):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < 1e-13:
            break
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    return nu

print(f"\n[2/6] Aplicando Pourbaix2016 (K_B={K_B} m/s)")
nu = true_anomaly(df['time'].values, P_yr_orb*365.25, e_orb, T0_BJD)
rv_model = V0_ms + K_B * (e_orb*np.cos(omega_B) + np.cos(omega_B + nu)) + dVB
df['rv_resid'] = df['rv'] - rv_model
df['rv_resid_centered'] = df['rv_resid'] - df['rv_resid'].median()

print(f"      rv mediana cruda:     {df['rv'].median():.2f} m/s")
print(f"      rv_model mediana:     {pd.Series(rv_model).median():.2f} m/s")
print(f"      rv_resid mediana:     {df['rv_resid'].median():.2f} m/s")
print(f"      rv_resid_centered std: {df['rv_resid_centered'].std():.3f} m/s")

# 3. Validar std
print(f"\n[3/6] VALIDACION: std post-Pourbaix por sub-segmento")
print(f"      {'segment':<35} {'N':>5} {'std_resid_centered':>20}")
for seg in sorted(df['model_inst'].unique()):
    m = df['model_inst'] == seg
    sub = df.loc[m, 'rv_resid_centered']
    sub_centered = sub - sub.median()
    std_seg = sub_centered.std()
    print(f"      {seg:<35} {m.sum():>5} {std_seg:>17.3f} m/s")

# Centrar por sub-segmento (offset tipico de instrumento)
df['rv_resid_centered_byseg'] = df.groupby('model_inst')['rv_resid'].transform(
    lambda x: x - x.median())

print(f"\n      std combinado tras centrar por segmento: "
      f"{df['rv_resid_centered_byseg'].std():.3f} m/s")
print(f"      ESPRESSO19_HR11 nuestro (referencia): 2.03 m/s")

# 4. Filtro outliers MAD 5sigma
print(f"\n[4/6] FILTRO OUTLIERS MAD (5-sigma robusto)")
def robust_sigma(x):
    return 1.4826 * (x - x.median()).abs().median()

df['sigma_robust_seg'] = df.groupby('model_inst')['rv_resid_centered_byseg'].transform(
    robust_sigma)
df['is_outlier'] = df['rv_resid_centered_byseg'].abs() > 5 * df['sigma_robust_seg']

print(f"      Outliers detectados: {df['is_outlier'].sum()}/{len(df)} "
      f"({df['is_outlier'].sum()/len(df)*100:.1f}%)")
for seg in sorted(df['model_inst'].unique()):
    m = (df['model_inst']==seg) & df['is_outlier']
    print(f"        {seg}: {m.sum()}")

df_clean = df[~df['is_outlier']].copy()
print(f"      N tras filtro: {len(df_clean)}")
for seg in sorted(df_clean['model_inst'].unique()):
    m = df_clean['model_inst']==seg
    sub_centered = df_clean.loc[m, 'rv_resid_centered_byseg']
    print(f"        {seg}: N={m.sum()}, "
          f"std_clean={sub_centered.std():.3f} m/s")

# 5. Apilar a nightly (mediana por noche)
print(f"\n[5/6] APILANDO A NIGHTLY (mediana por noche por segmento)")
df_clean['bjd_int'] = df_clean['time'].astype(int)

agg = df_clean.groupby(['bjd_int', 'model_inst']).agg(
    bjd=('time', 'mean'),
    rv=('rv', 'mean'),
    rv_resid=('rv_resid_centered_byseg', 'mean'),
    rv_err=('rv_err', lambda x: np.sqrt(np.sum(x**2))/len(x)),  # combinacion en cuadratura
    n_exp=('time', 'count'),
    ccf_fwhm=('ccf_fwhm', 'mean'),
    ccf_bispan=('ccf_bispan', 'mean'),
    ccf_contrast=('ccf_contrast', 'mean'),
).reset_index()

agg = agg.rename(columns={'model_inst': 'sub_inst'})
agg = agg.sort_values('bjd').reset_index(drop=True)

print(f"      Noches totales: {len(agg)}")
for seg in sorted(agg['sub_inst'].unique()):
    m = agg['sub_inst']==seg
    print(f"        {seg}: N_nights={m.sum()}, "
          f"std_resid={agg.loc[m,'rv_resid'].std():.3f} m/s")

# Guardar
agg.to_csv(f"{OUTDIR}/espresso_external_nightly.csv", index=False)
df_clean.to_csv(f"{OUTDIR}/espresso_external_post_pourbaix.csv", index=False)
print(f"      [OK] {OUTDIR}/espresso_external_nightly.csv")

# 6. Comparacion con ESPRESSO19_HR11 nuestro
print(f"\n[6/6] CROSS-CHECK CON ESPRESSO19_HR11 EXISTENTE")

master_path = f"{OUTDIR}/master_v2_nightly_residuals_clean.csv"
if os.path.exists(master_path):
    master = pd.read_csv(master_path)
    esp19 = master[master['sub_inst']=='ESPRESSO19_HR11'].copy()
    print(f"      Master_v2 ESPRESSO19_HR11: N={len(esp19)}")
    print(f"        BJD: {esp19['bjd'].min():.0f} - {esp19['bjd'].max():.0f}")
    print(f"        std_resid: {esp19['rv_resid'].std():.3f} m/s")
    
    # Cross-match por noche
    esp19['bjd_int'] = esp19['bjd'].astype(int)
    overlap = pd.merge(
        esp19[['bjd_int', 'rv_resid', 'rv_err']].rename(
            columns={'rv_resid': 'rv_resid_master', 'rv_err': 'rv_err_master'}),
        agg[['bjd_int', 'rv_resid', 'rv_err', 'sub_inst']].rename(
            columns={'rv_resid': 'rv_resid_external', 'rv_err': 'rv_err_external'}),
        on='bjd_int', how='inner')
    
    print(f"\n      Noches solapadas (en ambos datasets): {len(overlap)}")
    if len(overlap) > 0:
        diff = overlap['rv_resid_master'] - overlap['rv_resid_external']
        print(f"        Diferencia rv_resid (master - external):")
        print(f"          mediana: {diff.median():+.3f} m/s")
        print(f"          std:     {diff.std():.3f} m/s")
        print(f"          rango:   [{diff.min():+.2f}, {diff.max():+.2f}] m/s")
        if abs(diff.median()) > 1.0:
            print(f"        ⚠ OFFSET SISTEMATICO entre pipelines")
        else:
            print(f"        ✓ Pipelines consistentes en noches solapadas")
else:
    print(f"      [WARN] master_v2 no encontrado: {master_path}")

# ============================================================
# Plot diagnostico
# ============================================================
print(f"\n[PLOT] Generando diagnostico...")

fig, axes = plt.subplots(3, 1, figsize=(14, 10))

# Panel 1: serie temporal RV cruda
ax = axes[0]
for seg in sorted(df['model_inst'].unique()):
    m = df['model_inst']==seg
    ax.plot(df.loc[m,'time'], df.loc[m,'rv'], '.', ms=2, alpha=0.4,
            label=f'{seg} (N={m.sum()})')
ax.set_ylabel('rv crudo (m/s)')
ax.set_title('Serie temporal RV crudo (con binario)')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Panel 2: residuos post-Pourbaix
ax = axes[1]
for seg in sorted(df_clean['sub_inst' if 'sub_inst' in df_clean.columns else 'model_inst'].unique()):
    seg_col = 'model_inst' if 'model_inst' in df_clean.columns else 'sub_inst'
    m = df_clean[seg_col]==seg
    ax.errorbar(df_clean.loc[m,'time'],
                df_clean.loc[m,'rv_resid_centered_byseg'],
                yerr=df_clean.loc[m,'rv_err'],
                fmt='.', ms=2, alpha=0.4, label=f'{seg}')
ax.axhline(0, color='black', lw=0.5)
ax.set_ylabel('rv_resid centrado (m/s)')
ax.set_title(f'Residuos post-Pourbaix2016 + filtro outliers '
             f'(std={df_clean["rv_resid_centered_byseg"].std():.2f} m/s)')
ax.legend(fontsize=8)
ax.set_ylim(-15, 15)
ax.grid(alpha=0.3)

# Panel 3: nightly residuos
ax = axes[2]
for seg in sorted(agg['sub_inst'].unique()):
    m = agg['sub_inst']==seg
    ax.errorbar(agg.loc[m,'bjd'], agg.loc[m,'rv_resid'],
                yerr=agg.loc[m,'rv_err'],
                fmt='o', ms=4, label=f'{seg} N_nights={m.sum()}')
ax.axhline(0, color='black', lw=0.5)
ax.set_ylabel('rv_resid nightly (m/s)')
ax.set_xlabel('BJD')
ax.set_title('Residuos nightly post-Pourbaix (PARA INTEGRACION FUTURA)')
ax.legend(fontsize=8)
ax.set_ylim(-8, 8)
ax.grid(alpha=0.3)

fig.tight_layout()
out = f"{OUTDIR}/espresso_external_diagnostic.png"
fig.savefig(out, dpi=120)
print(f"      [OK] {out}")

print(f"\n{'='*72}")
print("RESUMEN FINAL")
print('='*72)
print(f"  Input:  1333 obs ESPRESSO 2023-2024 (1.4 yr baseline)")
print(f"  Output: {len(df_clean)} obs limpias, {len(agg)} noches")
print(f"  std post-Pourbaix nightly: {agg['rv_resid'].std():.3f} m/s")
print(f"  CCF indicators: ccf_fwhm, ccf_bispan, ccf_contrast (PRESERVADOS)")
print(f"\n  DECISION DE INTEGRACION:")
if df_clean['rv_resid_centered_byseg'].std() < 5.0:
    print(f"  ✓ std < 5 m/s → CANDIDATO A INTEGRACION en master_v3 manana")
else:
    print(f"  ✗ std > 5 m/s → INVESTIGAR antes de integrar")
print('='*72)
