"""
Fase 3 paso 1 - Resta del movimiento binario Pourbaix & Boffin 2016.

Usa los parametros publicados (Tabla 1) y el modelo:
  V_B(t) = V0 + K_B * (e*cos(omega_B) + cos(omega_B + nu(t))) + dVB

Calcula K_B analiticamente a partir de M_A, M_B, P, i, e.

Salida:
  data_interim/qc/master_v2_nightly_residuals.csv  (con columna rv_resid)
  data_interim/qc/binary_model_diagnostic.png
"""

import os
import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# Parametros Pourbaix & Boffin 2016 (Tabla 1, columna HARPS+ESO)
# ============================================================
P_yr        = 79.91          # periodo orbital
P_s         = P_yr * 365.25 * 86400.0
e           = 0.524           # excentricidad
i_deg       = 79.32           # inclinacion
i_rad       = np.deg2rad(i_deg)
omega_B_deg = 232.3           # argumento periastron de B
omega_B     = np.deg2rad(omega_B_deg)
T_jyear     = 1955.66         # paso por periastron (anio Juliano)
V0_kms      = -22.390         # velocidad sistemica
V0_ms       = V0_kms * 1000.0
dVB_ms      = 329.0           # offset gravitacional/convectivo de B
M_A         = 1.133           # masa A en M_sun
M_B         = 0.972           # masa B en M_sun
M_total     = M_A + M_B

# Convertir T (1955.66 anio juliano) a BJD
# astropy: jyear es 'Julian year', conversion exacta
T_obj = Time(T_jyear, format='jyear', scale='tdb')
T0_BJD = T_obj.jd
print(f"[INFO] T0 (paso periastron) = {T_jyear} jyear = BJD {T0_BJD:.4f}")

# ============================================================
# Calculo de K_B (semi-amplitud RV de B respecto al baricentro AB)
# ============================================================
# K_B = (2 pi G M_A^3 / (P (M_A+M_B)^2))^(1/3) * sin(i) / sqrt(1-e^2)
G = 6.67430e-11    # m^3 kg^-1 s^-2
M_sun = 1.98892e30 # kg

# K_B en m/s (formula estandar para SB1/SB2)
K_B = (2*np.pi*G * (M_A * M_sun)**3 / (P_s * ((M_A + M_B)*M_sun)**2))**(1/3) \
      * np.sin(i_rad) / np.sqrt(1 - e**2)

print(f"[INFO] K_B calculado = {K_B:.3f} m/s")
print(f"       (esperable ~6500-7500 m/s para alpha Cen B)")

# ============================================================
# Cargar master nightly
# ============================================================
INPUT  = "data_interim/qc/master_v2_nightly.csv"
OUTPUT = "data_interim/qc/master_v2_nightly_residuals.csv"
DIAG   = "data_interim/qc/binary_model_diagnostic.png"

df = pd.read_csv(INPUT)
print(f"\n[INFO] Cargado {INPUT}: {len(df)} epocas, {df['sub_inst'].nunique()} sub-instrumentos")

bjd = df['bjd'].values
n = len(bjd)

# ============================================================
# Modelo orbital
# ============================================================
def true_anomaly(t, P_days, e, T0_jd):
    """Calcula la anomalia verdadera nu(t) resolviendo Kepler."""
    M = 2*np.pi * (t - T0_jd) / P_days
    M = np.mod(M, 2*np.pi)
    
    # Resolver Kepler M = E - e*sin(E) por Newton
    E = M.copy()
    for _ in range(100):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < 1e-12:
            break
    
    # Anomalia verdadera
    nu = 2 * np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    return nu

P_days = P_yr * 365.25
nu = true_anomaly(bjd, P_days, e, T0_BJD)

# Modelo de RV de B en m/s (segun Pourbaix Eq 2)
rv_binary = V0_ms + K_B * (e*np.cos(omega_B) + np.cos(omega_B + nu)) + dVB_ms

df['rv_binary_model'] = rv_binary
df['rv_resid'] = df['rv'] - rv_binary

# ============================================================
# Diagnostico por sub-instrumento
# ============================================================
print("\n" + "="*70)
print("DIAGNOSTICO POR SUB-INSTRUMENTO")
print("="*70)
for sub in df['sub_inst'].unique():
    m = df['sub_inst'] == sub
    rv_med = df.loc[m, 'rv'].median()
    mod_med = df.loc[m, 'rv_binary_model'].median()
    resid_med = df.loc[m, 'rv_resid'].median()
    resid_std = df.loc[m, 'rv_resid'].std()
    print(f"  {sub:25s}  N={m.sum():4d}  "
          f"rv_med={rv_med:9.1f}  model_med={mod_med:9.1f}  "
          f"resid_med={resid_med:9.1f}  resid_std={resid_std:6.1f}")

# Guardar
df.to_csv(OUTPUT, index=False)
print(f"\n[OK] {OUTPUT}")

# ============================================================
# Plot diagnostico (3 paneles: data + modelo, residuos por inst, residuos zoom)
# ============================================================
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

fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=False)

# Panel 1: data + modelo
for sub in ORDER:
    m = df['sub_inst'] == sub
    if m.sum() == 0:
        continue
    axes[0].errorbar(df.loc[m, 'bjd'], df.loc[m, 'rv'], yerr=df.loc[m, 'rv_err'],
                     fmt='o', ms=3, alpha=0.6, color=COLORS[sub], label=sub)
# Modelo continuo
t_dense = np.linspace(bjd.min()-100, bjd.max()+100, 5000)
nu_dense = true_anomaly(t_dense, P_days, e, T0_BJD)
rv_model_dense = V0_ms + K_B * (e*np.cos(omega_B) + np.cos(omega_B + nu_dense)) + dVB_ms
axes[0].plot(t_dense, rv_model_dense, 'k-', lw=1, alpha=0.7,
             label='Pourbaix2016 model')
axes[0].set_ylabel('RV (m/s)')
axes[0].set_title('Datos vs modelo binario Pourbaix & Boffin 2016')
axes[0].legend(fontsize=8, loc='upper right', ncol=3)
axes[0].grid(alpha=0.3)

# Panel 2: residuos sin centrar
for sub in ORDER:
    m = df['sub_inst'] == sub
    if m.sum() == 0:
        continue
    axes[1].errorbar(df.loc[m, 'bjd'], df.loc[m, 'rv_resid'], yerr=df.loc[m, 'rv_err'],
                     fmt='o', ms=3, alpha=0.6, color=COLORS[sub], label=sub)
axes[1].axhline(0, color='k', lw=0.5)
axes[1].set_ylabel('rv_resid (m/s)')
axes[1].set_title('Residuos despues de restar Pourbaix2016 (cada sub-inst tiene su offset)')
axes[1].legend(fontsize=8, loc='upper right', ncol=3)
axes[1].grid(alpha=0.3)

# Panel 3: residuos centrados por sub-inst (mediana 0)
for sub in ORDER:
    m = df['sub_inst'] == sub
    if m.sum() == 0:
        continue
    rv_centered = df.loc[m, 'rv_resid'] - df.loc[m, 'rv_resid'].median()
    axes[2].errorbar(df.loc[m, 'bjd'], rv_centered, yerr=df.loc[m, 'rv_err'],
                     fmt='o', ms=3, alpha=0.6, color=COLORS[sub], label=sub)
axes[2].axhline(0, color='k', lw=0.5)
axes[2].set_ylabel('rv_resid centrado (m/s)')
axes[2].set_xlabel('BJD')
axes[2].set_title('Residuos centrados (rv_resid - mediana_sub_inst). ESTO es lo que va al periodograma.')
axes[2].set_ylim(-100, 100)  # zoom para ver los m/s
axes[2].legend(fontsize=8, loc='upper right', ncol=3)
axes[2].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(DIAG, dpi=120)
print(f"[OK] {DIAG}")

# ============================================================
# Estadistica final
# ============================================================
print("\n" + "="*70)
print("RESIDUOS GLOBALES TRAS RESTAR BINARIA")
print("="*70)
for sub in ORDER:
    m = df['sub_inst'] == sub
    if m.sum() == 0:
        continue
    centered = df.loc[m, 'rv_resid'] - df.loc[m, 'rv_resid'].median()
    print(f"  {sub:25s}  std centrado = {centered.std():7.2f} m/s")
