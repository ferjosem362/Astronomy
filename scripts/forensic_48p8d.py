"""
Test forense del pico a ~48.8 d.

Cinco bloques:
  C1) Aliases de cadencia: relacion con 1/Prot ± n/yr
  C2) Window function en banda 48-50 d
  C3) Pico por sub-instrumento individualmente
  C4) Test por mitades temporales (2004-2014 vs 2014-2024)
  C5) Comparacion con Dumusque.txt original (cross-check independiente)

Veredicto final: REAL / ALIAS / INDETERMINADO
"""

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"

P_TARGET = 48.80
P_LO_REL = 0.97
P_HI_REL = 1.03
P_BAND_LO = P_TARGET * P_LO_REL
P_BAND_HI = P_TARGET * P_HI_REL

# Aliases candidatos a probar
P_ROT = 38.65
P_YR = 365.25

print(f"\n{'='*72}")
print(f"FORENSIC TEST: pico a {P_TARGET} d (banda {P_BAND_LO:.2f}-{P_BAND_HI:.2f})")
print('='*72)

# ============================================================
# C1) Test aliases conocidos
# ============================================================
print(f"\n[C1] ALIASES CONOCIDOS")
print(f"     Probando 1/P_target = 1/{P_TARGET} = {1/P_TARGET:.6f} c/d")
print(f"     vs combinaciones de 1/Prot={1/P_ROT:.6f} y 1/yr={1/P_YR:.6f}")
print()

f_target = 1/P_TARGET
f_rot = 1/P_ROT
f_yr = 1/P_YR

best_alias = None
best_diff = np.inf
for n_rot in [-3,-2,-1,0,1,2,3]:
    for n_yr in [-4,-3,-2,-1,0,1,2,3,4]:
        for sign in [1, -1]:
            f_test = abs(n_rot * f_rot + n_yr * f_yr * sign)
            if f_test < 0.001:
                continue
            P_test = 1/f_test
            if 30 < P_test < 70:
                diff = abs(f_test - f_target) / f_target * 100  # % en frecuencia
                if diff < 2 and abs(n_rot) + abs(n_yr) <= 4:
                    print(f"     ¿{n_rot}/Prot {'+' if sign>0 else '-'} {n_yr}/yr? "
                          f"-> P={P_test:.3f} d (diff={diff:.2f}%)")
                    if diff < best_diff:
                        best_diff = diff
                        best_alias = (n_rot, n_yr, sign, P_test)

if best_alias:
    n_r, n_y, s, P_a = best_alias
    print(f"\n     MEJOR MATCH: P={P_a:.3f} d via {n_r}/Prot {'+' if s>0 else '-'} {n_y}/yr "
          f"(diff frec={best_diff:.2f}%)")
    if best_diff < 0.5:
        print(f"     ** ALIAS PROBABLE (match excelente) **")
    elif best_diff < 1.0:
        print(f"     * Alias plausible *")
else:
    print(f"     No se encontro alias simple obvio")

# ============================================================
# Cargar datos
# ============================================================
df = pd.read_csv(INPUT)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
df = df.sort_values('bjd').reset_index(drop=True)

t_all = df['bjd'].values
y_all = df['rv_resid_centered'].values
dy_all = np.maximum(df['rv_err'].values, 0.1)

# Frecuencias para todo
P_min, P_max = 1.5, 5000.0
nf = 50000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0/freqs
in_band = (periods >= P_BAND_LO) & (periods <= P_BAND_HI)

# ============================================================
# C2) Window function
# ============================================================
print(f"\n[C2] WINDOW FUNCTION en banda {P_BAND_LO:.2f}-{P_BAND_HI:.2f} d")
ls_win = LombScargle(t_all, np.ones_like(t_all), 1.0,
                     normalization='standard', center_data=False, fit_mean=False)
power_win = ls_win.power(freqs)
pwr_win_max = power_win[in_band].max()
P_win_max = periods[in_band][np.argmax(power_win[in_band])]
print(f"     Window power max en banda: {pwr_win_max:.4f} a P={P_win_max:.3f} d")
if pwr_win_max > 0.10:
    print(f"     ** SAMPLING TIENE PICO FUERTE EN ESTA BANDA → ALIAS LIKELY **")
elif pwr_win_max > 0.03:
    print(f"     * Sampling tiene pico moderado en banda *")
else:
    print(f"     - Sampling limpio en banda")

# ============================================================
# C3) Pico por sub-instrumento
# ============================================================
print(f"\n[C3] PICO POR SUB-INSTRUMENTO")
print(f"     Buscando en banda {P_BAND_LO:.2f}-{P_BAND_HI:.2f} d")
print(f"\n     {'sub_inst':<22} {'N':>4}  {'P_peak':>9} {'power':>9} {'FAP':>10} {'sig':>6}")
print(f"     {'-'*60}")

sub_results = {}
SUBS_TO_CHECK = ['HARPS03', 'ES_zhao2017', 'CHIRON_zhao2017', 'HARPS15',
                 'ESPRESSO19_HR11']

for sub in SUBS_TO_CHECK:
    m = df['sub_inst'] == sub
    n_pts = m.sum()
    if n_pts < 10:
        print(f"     {sub:<22} {n_pts:>4}  (sin datos suficientes)")
        continue
    
    t_s = df.loc[m, 'bjd'].values
    y_s = df.loc[m, 'rv_resid_centered'].values
    dy_s = np.maximum(df.loc[m, 'rv_err'].values, 0.1)
    
    ls_s = LombScargle(t_s, y_s, dy_s, normalization='standard')
    pwr_s = ls_s.power(freqs)
    
    pwr_in_band = pwr_s[in_band]
    if len(pwr_in_band) == 0:
        continue
    i_peak = np.argmax(pwr_in_band)
    P_peak = periods[in_band][i_peak]
    pwr_peak = pwr_in_band[i_peak]
    fap = ls_s.false_alarm_probability(pwr_peak, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    
    print(f"     {sub:<22} {n_pts:>4}  {P_peak:>9.3f} {pwr_peak:>9.4f} {fap:>10.2e}  {sig:>4}")
    sub_results[sub] = (P_peak, pwr_peak, fap)

# Veredicto C3
n_sig = sum(1 for v in sub_results.values() if v[2] < 0.01)
print(f"\n     Sub-instrumentos con FAP<1%: {n_sig}/{len(sub_results)}")
if n_sig >= 3:
    print(f"     ** APARECE EN VARIOS INSTRUMENTOS → ROBUSTO **")
elif n_sig == 2:
    print(f"     * Aparece en 2 instrumentos (moderado) *")
elif n_sig == 1:
    print(f"     - Solo en 1 instrumento (debil)")
else:
    print(f"     - No aparece significativamente en ningun instrumento individual")

# ============================================================
# C4) Test por mitades temporales
# ============================================================
print(f"\n[C4] TEST POR MITADES TEMPORALES")
bjd_mid = (t_all.min() + t_all.max()) / 2
m_first = t_all < bjd_mid
m_second = ~m_first

for label, mask in [('Mitad temprana', m_first), ('Mitad tardia', m_second)]:
    t_h = t_all[mask]
    y_h = y_all[mask]
    dy_h = dy_all[mask]
    if len(t_h) < 30:
        continue
    
    ls_h = LombScargle(t_h, y_h, dy_h, normalization='standard')
    pwr_h = ls_h.power(freqs)
    pwr_h_band = pwr_h[in_band]
    i_peak = np.argmax(pwr_h_band)
    P_peak = periods[in_band][i_peak]
    pwr_peak = pwr_h_band[i_peak]
    fap = ls_h.false_alarm_probability(pwr_peak, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"     {label:<20} N={len(t_h):>4}, "
          f"BJD {t_h.min():.0f}-{t_h.max():.0f}, "
          f"P_peak={P_peak:.3f} d, FAP={fap:.2e} {sig}")

# ============================================================
# C5) Comparacion con Dumusque.txt (test independiente)
# ============================================================
print(f"\n[C5] CROSS-CHECK CON DUMUSQUE.TXT (HARPS pre-fibra independiente)")

DUM_INPUT = "data_interim/qc/harps_pre_master_dumusque.csv"
import os
if os.path.exists(DUM_INPUT):
    df_dum = pd.read_csv(DUM_INPUT)
    print(f"     N={len(df_dum)}, baseline={(df_dum['bjd'].max()-df_dum['bjd'].min())/365.25:.2f} yr")
    
    # No tiene residuos post-Pourbaix; los calculamos rapido
    # Modelo Pourbaix simplificado (mismo que en subtract_binary_pourbaix2016.py)
    # Usar parametros directamente
    P_yr_orb = 79.91
    P_s_orb = P_yr_orb * 365.25 * 86400.0
    e_orb = 0.524
    i_rad = np.deg2rad(79.32)
    omega_B = np.deg2rad(232.3)
    T0_BJD = 2435349.815
    V0_ms = -22390.0
    dVB = 329.0
    K_B = 5503.97
    
    def true_anom(t, P_d, e, T0):
        M = 2*np.pi*(t - T0)/P_d % (2*np.pi)
        E = M.copy()
        for _ in range(100):
            dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
            E -= dE
            if np.max(np.abs(dE)) < 1e-12: break
        nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
        return nu
    
    nu = true_anom(df_dum['bjd'].values, P_yr_orb*365.25, e_orb, T0_BJD)
    rv_model = V0_ms + K_B * (e_orb*np.cos(omega_B) + np.cos(omega_B + nu)) + dVB
    df_dum['rv_resid_centered'] = (df_dum['rv'] - rv_model) - (df_dum['rv'] - rv_model).median()
    
    # Filtrar outliers MAD
    mad = (df_dum['rv_resid_centered'] - df_dum['rv_resid_centered'].median()).abs().median()
    sigma_r = 1.4826 * mad
    keep = df_dum['rv_resid_centered'].abs() < 5 * sigma_r
    df_dum_clean = df_dum[keep]
    print(f"     Tras filtro outliers: N={len(df_dum_clean)} (descartados {keep.sum()-len(df_dum_clean)})")
    print(f"     std residuos limpio: {df_dum_clean['rv_resid_centered'].std():.2f} m/s")
    
    # GLS
    t_d = df_dum_clean['bjd'].values
    y_d = df_dum_clean['rv_resid_centered'].values
    dy_d = np.maximum(df_dum_clean['rv_err'].values, 0.1)
    
    ls_d = LombScargle(t_d, y_d, dy_d, normalization='standard')
    pwr_d = ls_d.power(freqs)
    pwr_d_band = pwr_d[in_band]
    i_peak = np.argmax(pwr_d_band)
    P_peak = periods[in_band][i_peak]
    pwr_peak = pwr_d_band[i_peak]
    fap = ls_d.false_alarm_probability(pwr_peak, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"     Pico en banda 48.8 d (Dumusque only): P={P_peak:.3f} d, FAP={fap:.2e} {sig}")

# ============================================================
# VEREDICTO FINAL
# ============================================================
print(f"\n{'='*72}")
print(f"VEREDICTO FINAL para pico a {P_TARGET} d:")
print('='*72)

evidence = []
if best_alias and best_diff < 0.5:
    evidence.append(f"  ALIAS: P={best_alias[3]:.3f} es {best_alias[0]}/Prot ± {best_alias[1]}/yr (diff {best_diff:.2f}%)")
if pwr_win_max > 0.10:
    evidence.append(f"  ALIAS: window function tiene pico {pwr_win_max:.3f} en banda")
if n_sig >= 3:
    evidence.append(f"  REAL: aparece en {n_sig} sub-instrumentos con FAP<1%")
elif n_sig <= 1:
    evidence.append(f"  ALIAS: solo {n_sig} sub-inst con FAP<1%, no robusto")

print("\nEvidencias:")
for e in evidence:
    print(e)

print()
