"""
Capa 1 — Ingesta del Dumusque.txt como master HARPS_pre.

Convierte el .txt tab-separado del paper Nature 2012 (459 epocas HARPS de
Alpha Cen B 2008-2011) en un CSV limpio con columnas estandarizadas para
nuestro pipeline.

Salida: data_interim/qc/harps_pre_master_dumusque.csv
"""

import os
import sys
import numpy as np
import pandas as pd

INPUT_TXT  = "data_raw/Dumusque.txt"
OUTPUT_CSV = "data_interim/qc/harps_pre_master_dumusque.csv"
DIAG_PNG   = "data_interim/qc/harps_pre_dumusque_diagnostic.png"

if not os.path.exists(INPUT_TXT):
    sys.exit(f"[ERROR] No encuentro {INPUT_TXT}")

# Cabeceras segun documentacion del supplement Nature 2012
COLNAMES = [
    'jdb', 'rv_kms', 'rv_err_kms', 'rv_err_realistic_kms',
    'bis_kms', 'fwhm_kms', 'photon_noise_kms',
    'logrhk', 'logrhk_err', 'logrhk_lowpass',
    'alpha_hours', 'delta_deg', 'seeing'
]

print(f"[INFO] Leyendo {INPUT_TXT}")
print(f"[INFO] Asumiendo {len(COLNAMES)} columnas: {COLNAMES}")

# El fichero tiene tabs como separador, pero hemos visto en linea 3
# que algunos valores de logrhk pueden venir pegados (sin tab).
# Estrategia: leer con sep='\t', detectar filas con menos columnas de las
# esperadas, y reparar a mano si se puede.

# Saltamos las 2 primeras lineas (cabecera + separador de guiones)
df_raw = pd.read_csv(INPUT_TXT, sep='\t', skiprows=2, header=None,
                     engine='python', on_bad_lines='warn')

print(f"[INFO] Filas leidas: {len(df_raw)}")
print(f"[INFO] Columnas leidas: {df_raw.shape[1]} (esperabamos {len(COLNAMES)})")

# Diagnostico de columnas
if df_raw.shape[1] != len(COLNAMES):
    print("[WARN] Numero de columnas no coincide. Investigando filas problematicas...")
    # Buscamos las que tienen NaN en posiciones inesperadas
    for i in range(min(5, len(df_raw))):
        row = df_raw.iloc[i].tolist()
        n_nan = sum(1 for x in row if pd.isna(x))
        print(f"  Fila {i}: {len(row)} valores, {n_nan} NaN")

    # Si tenemos UNA columna menos, es el caso del logrhk pegado.
    # Lo intentamos reparar.
    if df_raw.shape[1] == len(COLNAMES) - 1:
        print("[FIX] Falta 1 columna. Probable: logrhk + logrhk_err pegados.")
        # Renombramos provisional
        cols_provisional = COLNAMES[:7] + ['logrhk_pegado'] + COLNAMES[9:]
        df_raw.columns = cols_provisional
        # Separar logrhk_pegado en dos
        # Formato esperado: '-4.99200.0005' -> -4.9920 y 0.0005
        def split_logrhk(v):
            if pd.isna(v):
                return (np.nan, np.nan)
            s = str(v)
            # logrhk siempre empieza con '-' y tiene 6 caracteres "-4.XXXX"
            if s.startswith('-') and len(s) >= 7:
                return (float(s[:7]), float(s[7:]))
            return (np.nan, np.nan)
        pairs = df_raw['logrhk_pegado'].apply(split_logrhk)
        df_raw['logrhk'] = [p[0] for p in pairs]
        df_raw['logrhk_err'] = [p[1] for p in pairs]
        df_raw = df_raw.drop(columns=['logrhk_pegado'])
        # Reordenar
        df_raw = df_raw[COLNAMES]
        print(f"[FIX] Reparado. Columnas finales: {list(df_raw.columns)}")
else:
    df_raw.columns = COLNAMES

df = df_raw.copy()

# ---------------------------------------------------------------
# Conversiones a unidades estandar de nuestro pipeline
# ---------------------------------------------------------------
print("\n[INFO] Aplicando conversiones de unidades...")

# RV km/s -> m/s
df['rv'] = df['rv_kms'] * 1000.0
df['rv_err'] = df['rv_err_realistic_kms'] * 1000.0   # usamos la realista
df['rv_err_formal'] = df['rv_err_kms'] * 1000.0
df['bis'] = df['bis_kms'] * 1000.0
df['fwhm'] = df['fwhm_kms'] * 1000.0
df['photon_noise'] = df['photon_noise_kms'] * 1000.0

# JDB -> BJD (jdb es JD-2400000 segun convencion HARPS)
df['bjd'] = df['jdb'] + 2400000.0

# Etiqueta de instrumento
df['instrument'] = 'HARPS_pre'
df['inst_split'] = 'HARPS_pre'

# ---------------------------------------------------------------
# QC
# ---------------------------------------------------------------
print(f"\n[STAT] N total: {len(df)}")
print(f"[STAT] Rango BJD: {df['bjd'].min():.1f} -> {df['bjd'].max():.1f}")
print(f"[STAT] Baseline: {df['bjd'].max() - df['bjd'].min():.1f} dias = "
      f"{(df['bjd'].max() - df['bjd'].min())/365.25:.2f} anos")
print(f"[STAT] RV mediana: {df['rv'].median():.3f} m/s")
print(f"[STAT] RV std:     {df['rv'].std():.3f} m/s")
print(f"[STAT] eRV mediana (realistic): {df['rv_err'].median():.3f} m/s")
print(f"[STAT] eRV mediana (formal):    {df['rv_err_formal'].median():.3f} m/s")
print(f"[STAT] FWHM mediana: {df['fwhm'].median():.1f} m/s")
print(f"[STAT] BIS mediana:  {df['bis'].median():.3f} m/s")
print(f"[STAT] log(R'HK) mediana: {df['logrhk'].median():.4f}")
print(f"[STAT] Seeing mediana: {df['seeing'].median():.2f} arcsec")
print(f"[STAT] Seeing rango:   {df['seeing'].min():.2f} -> {df['seeing'].max():.2f}")

# Aviso si hay NaNs
n_nan = df[['bjd','rv','rv_err','bis','fwhm','logrhk','seeing']].isna().sum()
if n_nan.sum() > 0:
    print(f"\n[WARN] NaNs detectados:")
    print(n_nan[n_nan > 0])

# ---------------------------------------------------------------
# Reordenar columnas a un esquema canonico para el resto del pipeline
# ---------------------------------------------------------------
out_cols = ['bjd', 'rv', 'rv_err', 'rv_err_formal',
            'bis', 'fwhm', 'photon_noise',
            'logrhk', 'logrhk_err', 'logrhk_lowpass',
            'seeing', 'alpha_hours', 'delta_deg',
            'instrument', 'inst_split']
df_out = df[out_cols].copy()

os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
df_out.to_csv(OUTPUT_CSV, index=False)
print(f"\n[OK] Escrito {OUTPUT_CSV}")
print(f"[OK] Columnas finales: {list(df_out.columns)}")

# ---------------------------------------------------------------
# Diagnostico visual
# ---------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

axes[0].errorbar(df_out['bjd'], df_out['rv'], yerr=df_out['rv_err'],
                 fmt='.', ms=3, alpha=0.6, color='tab:blue')
axes[0].set_ylabel('RV (m/s)')
axes[0].set_title(f"Dumusque.txt ingest: {len(df_out)} epocas HARPS_pre Alpha Cen B")
axes[0].grid(alpha=0.3)

axes[1].plot(df_out['bjd'], df_out['logrhk'], '.', ms=3, alpha=0.6, color='tab:green')
axes[1].plot(df_out['bjd'], df_out['logrhk_lowpass'], '-', lw=1, color='red',
             label='Low-pass filter (mag cycle)')
axes[1].set_ylabel("log(R'HK)")
axes[1].legend()
axes[1].grid(alpha=0.3)

axes[2].plot(df_out['bjd'], df_out['fwhm'], '.', ms=3, alpha=0.6, color='tab:orange')
axes[2].set_ylabel('FWHM (m/s)')
axes[2].grid(alpha=0.3)

axes[3].plot(df_out['bjd'], df_out['seeing'], '.', ms=3, alpha=0.6, color='tab:purple')
axes[3].axhline(2.20, color='gray', ls=':', lw=1, label='Seeing cuts Dumusque (2008-2011)')
axes[3].axhline(1.87, color='gray', ls=':', lw=1)
axes[3].axhline(1.63, color='gray', ls=':', lw=1)
axes[3].axhline(1.43, color='gray', ls=':', lw=1)
axes[3].set_ylabel('Seeing (arcsec)')
axes[3].set_xlabel('BJD')
axes[3].legend(fontsize=8)
axes[3].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(DIAG_PNG, dpi=120)
print(f"[OK] Figura {DIAG_PNG}")
