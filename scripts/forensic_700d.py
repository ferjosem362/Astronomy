"""
Test forense del pico a ~700 d.

Cinco bloques (adaptados a periodo largo):
  C1) Aliases candidatos (1yr, 2yr, ciclo magnetico, P_orb)
  C2) Window function en banda 650-750 d
  C3) Pico por sub-instrumento individual (excluyendo ESPRESSO19 por baseline)
  C4) Test por TERCIOS temporales (3 ciclos cada uno)
  C5) Cross-check informativo con Dumusque (limitado por baseline)
  C6) BONUS: refinar P_peak preciso con ajuste sinusoidal
"""

import os
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"

P_TARGET = 700.0
P_BAND_LO, P_BAND_HI = 600.0, 800.0

# Periodos relacionados a probar
P_YR     = 365.25
P_2YR    = 730.5
P_MAG_CYCLE_B = 7.8 * 365.25       # ciclo magnetico Alpha Cen B (~7.8 yr)
P_ORB_AB = 79.91 * 365.25          # binaria
P_ORB_AB_HARM = [P_ORB_AB/n for n in range(2, 50)]  # armonicos

print(f"\n{'='*72}")
print(f"FORENSIC TEST: pico a {P_TARGET} d (banda {P_BAND_LO}-{P_BAND_HI})")
print('='*72)

# ============================================================
# C1) Aliases candidatos para periodos largos
# ============================================================
print(f"\n[C1] ALIASES Y RELACIONES PARA PERIODOS LARGOS")
print(f"     (a) 1 yr = {P_YR:.2f} d")
print(f"     (b) 2 yr = {P_2YR:.2f} d")
print(f"     (c) Ciclo mag B (~7.8 yr) = {P_MAG_CYCLE_B:.0f} d")
print(f"     (d) Binaria AB = {P_ORB_AB:.0f} d")

# 700 d vs 2 yr: diferencia
diff_2yr = abs(P_TARGET - P_2YR) / P_2YR * 100
print(f"\n     P_target={P_TARGET} d vs 2*P_yr={P_2YR}: diff={diff_2yr:.1f}%")
if diff_2yr < 5:
    print(f"     ⚠ Cerca del alias bianual del sampling")

# Armonicos del binario que caen en banda
print(f"\n     Armonicos P_orb_AB / N en banda {P_BAND_LO}-{P_BAND_HI}:")
for n in range(2, 60):
    P_h = P_ORB_AB / n
    if P_BAND_LO < P_h < P_BAND_HI:
        diff = abs(P_h - P_TARGET) / P_TARGET * 100
        print(f"       N={n}: P={P_h:.1f} d  (diff con target={diff:.1f}%)")

# Cargar datos
df = pd.read_csv(INPUT)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
df = df.sort_values('bjd').reset_index(drop=True)

t_all = df['bjd'].values
y_all = df['rv_resid_centered'].values
dy_all = np.maximum(df['rv_err'].values, 0.1)

# Frecuencias
P_min, P_max = 1.5, 5000.0
nf = 50000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0/freqs
in_band = (periods >= P_BAND_LO) & (periods <= P_BAND_HI)

# ============================================================
# C2) Window function en banda
# ============================================================
print(f"\n[C2] WINDOW FUNCTION en banda {P_BAND_LO}-{P_BAND_HI} d")
ls_win = LombScargle(t_all, np.ones_like(t_all), 1.0,
                     normalization='standard', center_data=False, fit_mean=False)
power_win = ls_win.power(freqs)
pwr_win_max = power_win[in_band].max()
P_win_max = periods[in_band][np.argmax(power_win[in_band])]
print(f"     Window power max: {pwr_win_max:.4f} a P={P_win_max:.1f} d")
if pwr_win_max > 0.30:
    print(f"     ⚠ ALIAS LIKELY: window function tiene pico fuerte")
elif pwr_win_max > 0.10:
    print(f"     * Window moderado, posible alias")
else:
    print(f"     ✓ Window relativamente limpia en banda")

# ============================================================
# C3) Pico por sub-instrumento INDIVIDUAL
# ============================================================
print(f"\n[C3] PICO POR SUB-INSTRUMENTO (excluyendo ESPRESSO19 por baseline)")
print(f"\n     {'sub_inst':<22} {'N':>4} {'baseline_yr':>12} "
      f"{'P_peak':>9} {'power':>9} {'FAP':>10} {'sig':>6} {'ciclos':>8}")
print(f"     {'-'*92}")

# Solo sub-inst con baseline >= 1.5 ciclos del periodo target
SUBS_TO_CHECK = ['HARPS03', 'ES_zhao2017', 'CHIRON_zhao2017', 'HARPS15']

sub_results = {}
for sub in SUBS_TO_CHECK:
    m = df['sub_inst'] == sub
    n_pts = m.sum()
    if n_pts < 10:
        print(f"     {sub:<22} {n_pts:>4}  insuficiente")
        continue
    
    t_s = df.loc[m, 'bjd'].values
    baseline = (t_s.max() - t_s.min()) / 365.25
    n_cycles = (t_s.max() - t_s.min()) / P_TARGET
    
    if n_cycles < 1.0:
        print(f"     {sub:<22} {n_pts:>4} {baseline:>12.2f} "
              f"  baseline insuficiente ({n_cycles:.1f} ciclos)")
        continue
    
    y_s = df.loc[m, 'rv_resid_centered'].values
    dy_s = np.maximum(df.loc[m, 'rv_err'].values, 0.1)
    
    ls_s = LombScargle(t_s, y_s, dy_s, normalization='standard')
    pwr_s = ls_s.power(freqs)
    pwr_in_band = pwr_s[in_band]
    i_peak = np.argmax(pwr_in_band)
    P_peak = periods[in_band][i_peak]
    pwr_peak = pwr_in_band[i_peak]
    fap = ls_s.false_alarm_probability(pwr_peak, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    
    print(f"     {sub:<22} {n_pts:>4} {baseline:>12.2f} "
          f"{P_peak:>9.1f} {pwr_peak:>9.4f} {fap:>10.2e}  {sig:>4}  {n_cycles:>5.1f}")
    sub_results[sub] = (P_peak, pwr_peak, fap, n_cycles)

# Dejamos ESPRESSO19 fuera del veredicto pero lo reportamos
print(f"\n     [Para info, NO contado en veredicto:]")
m = df['sub_inst'] == 'ESPRESSO19_HR11'
if m.sum() > 0:
    t_s = df.loc[m, 'bjd'].values
    baseline = (t_s.max() - t_s.min()) / 365.25
    n_cycles = (t_s.max() - t_s.min()) / P_TARGET
    y_s = df.loc[m, 'rv_resid_centered'].values
    dy_s = np.maximum(df.loc[m, 'rv_err'].values, 0.1)
    ls_s = LombScargle(t_s, y_s, dy_s, normalization='standard')
    pwr_s = ls_s.power(freqs)
    i_peak = np.argmax(pwr_s[in_band])
    P_peak = periods[in_band][i_peak]
    pwr_peak = pwr_s[in_band][i_peak]
    fap = ls_s.false_alarm_probability(pwr_peak, method='baluev')
    print(f"     ESPRESSO19_HR11 (baseline {baseline:.2f}yr={n_cycles:.2f} ciclos): "
          f"P={P_peak:.1f}, FAP={fap:.2e}  [TRAMPA: pendiente lineal]")

n_sig = sum(1 for v in sub_results.values() if v[2] < 0.05)
print(f"\n     Sub-instrumentos con FAP<5%: {n_sig}/{len(sub_results)}")

# ============================================================
# C4) Test por TERCIOS temporales
# ============================================================
print(f"\n[C4] TEST POR TERCIOS TEMPORALES")
bjd_min, bjd_max = t_all.min(), t_all.max()
range_full = bjd_max - bjd_min
edges = [bjd_min + i*range_full/3 for i in range(4)]

print(f"     Rango total: {bjd_min:.0f} - {bjd_max:.0f} ({range_full/365.25:.1f} yr)")
for i in range(3):
    mask = (t_all >= edges[i]) & (t_all <= edges[i+1])
    t_t = t_all[mask]; y_t = y_all[mask]; dy_t = dy_all[mask]
    if len(t_t) < 30:
        print(f"     Tercio {i+1} N={len(t_t)} insuficiente")
        continue
    
    n_cycles_t = (t_t.max() - t_t.min()) / P_TARGET
    ls_t = LombScargle(t_t, y_t, dy_t, normalization='standard')
    pwr_t = ls_t.power(freqs)
    pwr_t_band = pwr_t[in_band]
    i_peak = np.argmax(pwr_t_band)
    P_peak = periods[in_band][i_peak]
    pwr_peak = pwr_t_band[i_peak]
    fap = ls_t.false_alarm_probability(pwr_peak, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"     Tercio {i+1} ({edges[i]:.0f}-{edges[i+1]:.0f}, N={len(t_t)}, "
          f"{n_cycles_t:.1f} ciclos): P={P_peak:.1f} d, FAP={fap:.2e} {sig}")

# ============================================================
# C5) Cross-check Dumusque (limitado pero informativo)
# ============================================================
print(f"\n[C5] CROSS-CHECK CON DUMUSQUE.TXT (3.37 yr = {3.37*365.25/P_TARGET:.1f} ciclos de {P_TARGET}d)")

DUM_INPUT = "data_interim/qc/harps_pre_master_dumusque.csv"
if os.path.exists(DUM_INPUT):
    df_dum = pd.read_csv(DUM_INPUT)
    
    # Modelo Pourbaix
    P_yr_orb = 79.91
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
    rv_model = V0_ms + K_B*(e_orb*np.cos(omega_B) + np.cos(omega_B + nu)) + dVB
    df_dum['rv_resid_centered'] = (df_dum['rv'] - rv_model) - (df_dum['rv'] - rv_model).median()
    
    mad = (df_dum['rv_resid_centered'] - df_dum['rv_resid_centered'].median()).abs().median()
    sigma_r = 1.4826 * mad
    keep = df_dum['rv_resid_centered'].abs() < 5*sigma_r
    df_dum = df_dum[keep]
    
    t_d = df_dum['bjd'].values
    y_d = df_dum['rv_resid_centered'].values
    dy_d = np.maximum(df_dum['rv_err'].values, 0.1)
    
    ls_d = LombScargle(t_d, y_d, dy_d, normalization='standard')
    pwr_d = ls_d.power(freqs)
    i_peak = np.argmax(pwr_d[in_band])
    P_peak = periods[in_band][i_peak]
    pwr_peak = pwr_d[in_band][i_peak]
    fap = ls_d.false_alarm_probability(pwr_peak, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"     N={len(df_dum)}, baseline 3.37 yr (limitado para 700d)")
    print(f"     Pico mas alto en banda: P={P_peak:.1f} d, FAP={fap:.2e} {sig}")
    print(f"     Nota: con 1.7 ciclos no es test definitivo, pero informa direccion")

# ============================================================
# C6) Refinar P y K con ajuste sinusoidal
# ============================================================
print(f"\n[C6] REFINAMIENTO DE P, K POR AJUSTE SINUSOIDAL")
ls_full = LombScargle(t_all, y_all, dy_all, normalization='standard')
power_full = ls_full.power(freqs)
i_peak_glob = np.argmax(power_full * in_band)
P_peak_glob = periods[in_band][np.argmax(power_full[in_band])]
print(f"     P_peak GLS combinado: {P_peak_glob:.2f} d")

# Ajuste sinusoidal fino con curve_fit
def sinusoid(t, K, P, phi, offset):
    return K*np.sin(2*np.pi*t/P + phi) + offset

try:
    p0 = [3.0, P_peak_glob, 0.0, 0.0]
    bounds = ([0, P_peak_glob*0.95, -2*np.pi, -10],
              [50, P_peak_glob*1.05, 2*np.pi, 10])
    popt, pcov = curve_fit(sinusoid, t_all, y_all, p0=p0,
                           sigma=dy_all, absolute_sigma=True, bounds=bounds)
    K_fit, P_fit, phi_fit, off_fit = popt
    K_err, P_err, _, _ = np.sqrt(np.diag(pcov))
    
    # Calcular masa equivalente
    M_STAR = 0.972
    P_yr_calc = P_fit / 365.25
    M_p_earth = K_fit / 0.6395 * (M_STAR)**(2/3) * (P_yr_calc)**(1/3)
    
    print(f"     P_fit  = {P_fit:.2f} ± {P_err:.2f} d")
    print(f"     K_fit  = {K_fit:.3f} ± {K_err:.3f} m/s")
    print(f"     M·sini = {M_p_earth:.2f} M_Earth (asumiendo i=90, e=0)")
    
    # Calcular chi2 reducido
    y_model = sinusoid(t_all, *popt)
    resid_model = y_all - y_model
    chi2 = np.sum((resid_model/dy_all)**2)
    chi2_red = chi2 / (len(t_all) - 4)
    print(f"     chi2_red = {chi2_red:.2f}")
except Exception as e:
    print(f"     Error en ajuste: {e}")
    K_fit, P_fit = None, None

# ============================================================
# VEREDICTO
# ============================================================
print(f"\n{'='*72}")
print(f"VEREDICTO PRELIMINAR para pico a ~{P_TARGET} d:")
print('='*72)

evidence_pro = []
evidence_con = []

if pwr_win_max < 0.10:
    evidence_pro.append(f"  Window function limpia en banda ({pwr_win_max:.3f})")
elif pwr_win_max > 0.30:
    evidence_con.append(f"  Window con pico fuerte ({pwr_win_max:.3f})")

if n_sig >= 2:
    evidence_pro.append(f"  Aparece en {n_sig} sub-inst con FAP<5% (sin contar ESPRESSO19)")
else:
    evidence_con.append(f"  Solo en {n_sig} sub-inst con FAP<5%")

# El dato del HARPS03 individual con FAP<1% en GLS combinado original es positivo
# Pero hay que verificar consistencia entre tercios

print("\nA FAVOR (señal real):")
for e in evidence_pro:
    print(e)
print("\nEN CONTRA (artefacto):")
for e in evidence_con:
    print(e)

print(f"\nSe necesita interpretacion humana mirando tabla C3, C4, C5 y K_fit")
