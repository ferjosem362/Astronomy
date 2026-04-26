"""
PASO 1: Corregir offset de 39.6 m/s en las 52 noches externas
        -> master_v3_v2_nightly_residuals_clean.csv
PASO 2: Test forense del pico a ~232 d (similar al de 48.8 d y 700 d)
        Cinco bloques: aliases, window, sub-instrumentos, tercios, ajuste
"""

import os
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_V3   = "data_interim/qc/master_v3_nightly_residuals_clean.csv"
INPUT_EXT  = "data_interim/qc/espresso_external_nightly.csv"
OUTDIR     = "data_interim/qc"
OUTPUT_V3v2 = f"{OUTDIR}/master_v3_v2_nightly_residuals_clean.csv"

OFFSET = 39.634  # m/s, mediana medida en solapamiento ESPRESSO19_HR11 vs externo

P_TARGET = 232.26
P_LO = P_TARGET * 0.97
P_HI = P_TARGET * 1.03

# ============================================================
# PASO 1: corregir offset
# ============================================================
print("\n" + "="*72)
print(f"PASO 1: APLICAR OFFSET -{OFFSET} m/s A 52 NOCHES EXTERNAS")
print("="*72)

df_v3 = pd.read_csv(INPUT_V3)
print(f"\n  Master_v3 original: {len(df_v3)} noches")

# Sub-instrumentos del externo a corregir
EXT_SUBS = ['ESPRESSO_EXT_2023e', 'ESPRESSO_EXT_2023plus']
mask_ext = df_v3['sub_inst'].isin(EXT_SUBS)
n_ext = mask_ext.sum()
print(f"  Noches externas (a corregir): {n_ext}")

# Aplicar offset (sumamos +offset a las externas para alinearlas con v2 que estaba "ext + offset")
# Pero diff_v2_ext = -39.6 → v2 = ext - 39.6 → ext = v2 + 39.6
# Para alinear ext al sistema de v2, restamos -39.6 (i.e., sumamos -39.6 a las RV externas)
# Equivalente: rv_resid_ext_corrected = rv_resid_ext - 39.6
# Wait: la diferencia v2 - ext = -39.6 significa v2 está -39.6 debajo de ext
#       Para que ext se alinee con v2: ext_new = ext + (v2 - ext) = ext - 39.6
# Aplicamos -39.6 a rv_resid del externo
df_v3.loc[mask_ext, 'rv_resid'] = df_v3.loc[mask_ext, 'rv_resid'] - OFFSET

# Recalcular rv_resid_centered por sub_inst
df_v3['rv_resid_centered'] = df_v3.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())

# Verificar: ahora el solapamiento debería ser cercano a cero
df_v3_sorted = df_v3.sort_values('bjd').reset_index(drop=True)
df_v3_sorted.to_csv(OUTPUT_V3v2, index=False)
print(f"  [OK] {OUTPUT_V3v2}")

# Verificacion del offset corregido
print(f"\n  Verificacion: stat por sub_inst tras correccion")
print(f"  {'sub_inst':<28} {'N':>4}  {'std':>8}  {'med':>9}")
for sub in sorted(df_v3['sub_inst'].unique()):
    m = df_v3['sub_inst'] == sub
    sub_data = df_v3.loc[m, 'rv_resid']
    print(f"  {sub:<28} {m.sum():>4}  {sub_data.std():>7.2f}  {sub_data.median():>+8.2f}")

# ============================================================
# PASO 2: TEST FORENSE 232 d
# ============================================================
print("\n" + "="*72)
print(f"PASO 2: FORENSIC TEST pico ~{P_TARGET} d (banda {P_LO:.2f}-{P_HI:.2f})")
print("="*72)

# Cargar master_v3_v2 (con offset corregido)
df = pd.read_csv(OUTPUT_V3v2)
t_all = df['bjd'].values
y_all = df['rv_resid_centered'].values
dy_all = np.maximum(df['rv_err'].values, 0.1)

# Frecuencias
P_min, P_max = 1.5, 5000.0
nf = 50000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0/freqs
in_band = (periods >= P_LO) & (periods <= P_HI)

# ============================================================
# C1) Aliases conocidos
# ============================================================
print(f"\n[C1] ALIASES Y RELACIONES")

P_ROT = 38.65
P_YR  = 365.25

f_target = 1.0/P_TARGET
print(f"     1/P_target = {f_target:.6f} c/d")
print(f"     vs combinaciones de 1/Prot, 1/yr, 2/yr, ...")

# Probar múltiplos enteros de Prot directos
print(f"\n     Múltiplos de Prot=38.65 d:")
for n in range(2, 12):
    P_n = P_ROT * n
    if 200 < P_n < 280:
        diff_pct = abs(P_n - P_TARGET) / P_TARGET * 100
        sig = "** ARMONICO PROBABLE **" if diff_pct < 1 else ""
        print(f"       {n}*Prot = {P_n:.2f} d  (diff con target {diff_pct:.2f}%) {sig}")

# Combinación 1/Prot ± n/yr
print(f"\n     Combinaciones 1/Prot ± n/yr:")
best_alias = None
best_diff = np.inf
for n_rot in [1, 2, 3]:
    for n_yr in range(-5, 6):
        for sign in [1, -1]:
            f_test = abs(1.0/(P_ROT/n_rot) + n_yr * (1.0/P_YR) * sign)
            if f_test < 0.0001: continue
            P_test = 1.0/f_test
            if 200 < P_test < 280:
                diff_pct = abs(f_test - f_target) / f_target * 100
                if diff_pct < 5:
                    print(f"       {n_rot}/Prot {'+' if sign>0 else '-'} {n_yr}/yr -> "
                          f"P={P_test:.2f} d  (diff {diff_pct:.2f}%)")
                    if diff_pct < best_diff:
                        best_diff = diff_pct
                        best_alias = (n_rot, n_yr, sign, P_test)

# 1/yr armonicos
print(f"\n     Armónicos de 1/yr:")
for n in [1, 2, 3, 4, 5]:
    P_n = P_YR / n
    diff_pct = abs(P_n - P_TARGET) / P_TARGET * 100
    if diff_pct < 10:
        print(f"       1/{n}*yr = {P_n:.2f} d  (diff {diff_pct:.2f}%)")
    P_n2 = P_YR * n / (n-0.5) if n > 0.5 else 1e9
    if 200 < P_n2 < 280:
        diff_pct = abs(P_n2 - P_TARGET) / P_TARGET * 100
        print(f"       n*yr / (n-0.5) [n={n}] = {P_n2:.2f}  (diff {diff_pct:.2f}%)")

# ============================================================
# C2) Window function
# ============================================================
print(f"\n[C2] WINDOW FUNCTION en banda {P_LO:.2f}-{P_HI:.2f}")
ls_win = LombScargle(t_all, np.ones_like(t_all), 1.0,
                     normalization='standard', center_data=False, fit_mean=False)
power_win = ls_win.power(freqs)
pwr_win_max = power_win[in_band].max()
P_win_max = periods[in_band][np.argmax(power_win[in_band])]
print(f"     Window max en banda: {pwr_win_max:.4f} a P={P_win_max:.2f} d")
if pwr_win_max > 0.30:
    print(f"     ⚠ ALIAS LIKELY")
elif pwr_win_max > 0.10:
    print(f"     * Window moderado")
else:
    print(f"     ✓ Window limpia")

# ============================================================
# C3) Pico por sub-instrumento
# ============================================================
print(f"\n[C3] PICO POR SUB-INSTRUMENTO INDIVIDUAL")
print(f"     {'sub_inst':<28} {'N':>4}  {'baseline_yr':>11} {'P_peak':>8} "
      f"{'power':>8} {'FAP':>10} {'sig':>5} {'n_cyc':>6}")
print(f"     {'-'*92}")

SUBS = ['HARPS03', 'ES_zhao2017', 'CHIRON_zhao2017', 'HARPS15',
        'ESPRESSO19_HR11', 'ESPRESSO_EXT_2023plus', 'ESPRESSO_EXT_2023e']

sub_results = {}
for sub in SUBS:
    m = df['sub_inst'] == sub
    n_pts = m.sum()
    if n_pts < 8:
        continue
    
    t_s = df.loc[m, 'bjd'].values
    baseline = (t_s.max() - t_s.min()) / 365.25
    n_cyc = (t_s.max() - t_s.min()) / P_TARGET
    
    if n_cyc < 0.8:
        print(f"     {sub:<28} {n_pts:>4}  {baseline:>11.2f} "
              f"  baseline insuficiente ({n_cyc:.1f} ciclos)")
        continue
    
    y_s = df.loc[m, 'rv_resid_centered'].values
    dy_s = np.maximum(df.loc[m, 'rv_err'].values, 0.1)
    
    ls_s = LombScargle(t_s, y_s, dy_s, normalization='standard')
    pwr_s = ls_s.power(freqs)
    pwr_b = pwr_s[in_band]
    i_p = np.argmax(pwr_b)
    P_p = periods[in_band][i_p]
    pwr_p = pwr_b[i_p]
    fap = ls_s.false_alarm_probability(pwr_p, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"     {sub:<28} {n_pts:>4}  {baseline:>11.2f} {P_p:>8.2f} "
          f"{pwr_p:>8.4f} {fap:>10.2e}  {sig:>4}  {n_cyc:>5.1f}")
    sub_results[sub] = (P_p, pwr_p, fap, n_cyc)

n_sig = sum(1 for v in sub_results.values() if v[2] < 0.05)
print(f"\n     Sub-instrumentos con FAP<5%: {n_sig}/{len(sub_results)}")

# Verificar consistencia P_peak entre sub-inst
P_peaks = [v[0] for v in sub_results.values() if v[2] < 0.10]
if len(P_peaks) >= 2:
    P_relstd = np.std(P_peaks) / np.mean(P_peaks) * 100
    print(f"     P_peak std relativo entre sub-inst (FAP<10%): {P_relstd:.2f}%")
    if P_relstd < 1:
        print(f"     ✓ P_peak consistente al <1% — coherente con Kepleriano")
    elif P_relstd < 3:
        print(f"     * P_peak consistente al <3% — aceptable")
    else:
        print(f"     ⚠ P_peak migra >{P_relstd:.0f}% — sospechoso")

# ============================================================
# C4) Test por TERCIOS
# ============================================================
print(f"\n[C4] TEST POR TERCIOS TEMPORALES")
bjd_min, bjd_max = t_all.min(), t_all.max()
edges = [bjd_min + i*(bjd_max-bjd_min)/3 for i in range(4)]

for i in range(3):
    mask = (t_all >= edges[i]) & (t_all <= edges[i+1])
    t_t, y_t, dy_t = t_all[mask], y_all[mask], dy_all[mask]
    if len(t_t) < 30: continue
    n_cyc = (t_t.max() - t_t.min()) / P_TARGET
    
    ls_t = LombScargle(t_t, y_t, dy_t, normalization='standard')
    pwr_t = ls_t.power(freqs)
    pwr_b = pwr_t[in_band]
    i_p = np.argmax(pwr_b)
    P_p = periods[in_band][i_p]
    pwr_p = pwr_b[i_p]
    fap = ls_t.false_alarm_probability(pwr_p, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"     Tercio {i+1} ({edges[i]:.0f}-{edges[i+1]:.0f}, N={len(t_t)}, "
          f"{n_cyc:.1f} ciclos): P={P_p:.2f}, FAP={fap:.2e} {sig}")

# ============================================================
# C5) Ajuste sinusoidal y K, masa
# ============================================================
print(f"\n[C5] AJUSTE SINUSOIDAL")

def sinusoid(t, K, P, phi, offset):
    return K*np.sin(2*np.pi*t/P + phi) + offset

# Buscar pico exacto en banda
ls_full = LombScargle(t_all, y_all, dy_all, normalization='standard')
power_full = ls_full.power(freqs)
i_p = np.argmax(power_full * in_band)
P_init = periods[in_band][np.argmax(power_full[in_band])]
print(f"     P pico en banda: {P_init:.3f} d")

try:
    p0 = [2.0, P_init, 0.0, 0.0]
    bounds = ([0, P_init*0.95, -2*np.pi, -10],
              [50, P_init*1.05, 2*np.pi, 10])
    popt, pcov = curve_fit(sinusoid, t_all, y_all, p0=p0,
                           sigma=dy_all, absolute_sigma=True, bounds=bounds)
    K_fit, P_fit, _, _ = popt
    K_err, P_err, _, _ = np.sqrt(np.diag(pcov))
    
    M_STAR = 0.972
    P_yr_calc = P_fit / 365.25
    M_p_earth = K_fit / 0.6395 * (M_STAR)**(2/3) * (P_yr_calc)**(1/3)
    
    y_model = sinusoid(t_all, *popt)
    resid = y_all - y_model
    chi2 = np.sum((resid/dy_all)**2)
    chi2_red = chi2 / (len(t_all) - 4)
    
    print(f"     P_fit  = {P_fit:.3f} ± {P_err:.3f} d")
    print(f"     K_fit  = {K_fit:.3f} ± {K_err:.3f} m/s")
    print(f"     M·sini = {M_p_earth:.3f} M_Earth")
    print(f"     chi2_red = {chi2_red:.2f}")
    if chi2_red < 5:
        print(f"     ✓ Modelo sinusoidal AJUSTA bien")
    elif chi2_red < 50:
        print(f"     * chi2_red moderado, modelo aceptable pero hay residuo")
    else:
        print(f"     ⚠ chi2_red >>1, modelo NO describe los datos")
except Exception as e:
    print(f"     Error en ajuste: {e}")

# ============================================================
# C6) Cross-check Dumusque (independiente)
# ============================================================
print(f"\n[C6] CROSS-CHECK DUMUSQUE (HARPS pre-fibra independiente)")
DUM_INPUT = "data_interim/qc/harps_pre_master_dumusque.csv"
if os.path.exists(DUM_INPUT):
    df_dum = pd.read_csv(DUM_INPUT)
    
    # Aplicar Pourbaix
    P_yr_orb = 79.91; e_orb = 0.524
    omega_B = np.deg2rad(232.3); T0 = 2435349.815
    V0 = -22390.0; dVB = 329.0; K_B = 5503.97
    
    def true_anom(t, P_d, e, T0):
        M = 2*np.pi*(t - T0)/P_d % (2*np.pi)
        E = M.copy()
        for _ in range(100):
            dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
            E -= dE
            if np.max(np.abs(dE)) < 1e-12: break
        return 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    
    nu = true_anom(df_dum['bjd'].values, P_yr_orb*365.25, e_orb, T0)
    rv_model = V0 + K_B*(e_orb*np.cos(omega_B) + np.cos(omega_B + nu)) + dVB
    df_dum['rv_resid'] = (df_dum['rv'] - rv_model) - (df_dum['rv'] - rv_model).median()
    
    mad = (df_dum['rv_resid'] - df_dum['rv_resid'].median()).abs().median()
    sigma_r = 1.4826 * mad
    df_dum = df_dum[df_dum['rv_resid'].abs() < 5*sigma_r]
    
    n_cyc = (df_dum['bjd'].max() - df_dum['bjd'].min()) / P_TARGET
    print(f"     N={len(df_dum)}, baseline {(df_dum['bjd'].max()-df_dum['bjd'].min())/365.25:.2f} yr, "
          f"{n_cyc:.1f} ciclos de {P_TARGET}d")
    
    ls_d = LombScargle(df_dum['bjd'], df_dum['rv_resid'],
                       np.maximum(df_dum['rv_err'], 0.1),
                       normalization='standard')
    pwr_d = ls_d.power(freqs)
    pwr_d_b = pwr_d[in_band]
    i_p = np.argmax(pwr_d_b)
    P_p = periods[in_band][i_p]
    pwr_p = pwr_d_b[i_p]
    fap = ls_d.false_alarm_probability(pwr_p, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"     Pico en banda 232d: P={P_p:.2f}, FAP={fap:.2e} {sig}")

# ============================================================
# Veredicto
# ============================================================
print(f"\n{'='*72}")
print(f"VEREDICTO PRELIMINAR para pico {P_TARGET} d")
print('='*72)

# Plot diagnóstico
fig, axes = plt.subplots(2, 2, figsize=(15, 10))

# Zoom espectral
ax = axes[0,0]
P_view_lo, P_view_hi = P_TARGET*0.92, P_TARGET*1.08
mask_view = (periods >= P_view_lo) & (periods <= P_view_hi)
ax.plot(periods[mask_view], power_full[mask_view], 'k-', lw=1)
fap_th = ls_full.false_alarm_level([0.10,0.01,0.001], method='baluev')
for fl, fp in zip([10,1,0.1], fap_th):
    ax.axhline(fp, color='gray', ls='--', lw=0.6, alpha=0.6, label=f'FAP {fl}%')
ax.axvline(P_TARGET, color='green', ls=':', lw=2, label=f'P={P_TARGET}')
# Marcar 6*Prot
ax.axvline(6*P_ROT, color='red', ls=':', lw=1.5, alpha=0.7, label=f'6×Prot={6*P_ROT:.1f}')
ax.set_xlabel('Period (d)')
ax.set_ylabel('GLS power')
ax.set_title(f'Zoom espectral 232 d')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Power por sub-inst (bar plot)
ax = axes[0,1]
if sub_results:
    sub_names = list(sub_results.keys())
    sub_pows = [sub_results[s][1] for s in sub_names]
    sub_faps = [sub_results[s][2] for s in sub_names]
    colors_bar = ['green' if f<0.01 else ('orange' if f<0.10 else 'red') for f in sub_faps]
    sub_short = [s.replace('_zhao2017','_z').replace('_HR11','').replace('ESPRESSO_','ESP_') for s in sub_names]
    ax.bar(sub_short, sub_pows, color=colors_bar, edgecolor='black')
    ax.tick_params(axis='x', rotation=30)
    ax.set_ylabel('GLS power en banda')
    ax.set_title(f'Pico {P_TARGET}d por sub-inst (verde=FAP<1%)')
    ax.grid(alpha=0.3, axis='y')

# Phase folded
ax = axes[1,0]
phase = ((t_all - t_all[0]) / P_TARGET) % 1
ax.errorbar(phase, y_all, yerr=dy_all, fmt='.', ms=3, alpha=0.5)
ax.set_xlabel(f'Phase (P={P_TARGET}d)')
ax.set_ylabel('rv_resid (m/s)')
ax.set_title(f'Phase-folded a P={P_TARGET}d')
ax.set_ylim(-15, 15)
ax.grid(alpha=0.3)

# Window function en banda extendida
ax = axes[1,1]
P_view_lo, P_view_hi = 200, 280
mask_view = (periods >= P_view_lo) & (periods <= P_view_hi)
ax.plot(periods[mask_view], power_win[mask_view], 'darkred', lw=1)
ax.axvline(P_TARGET, color='green', ls=':', lw=2, label=f'P={P_TARGET}')
ax.set_xlabel('Period (d)')
ax.set_ylabel('Window power')
ax.set_title('Window function en banda 200-280 d')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

fig.tight_layout()
out_plot = f"{OUTDIR}/forensic_232d.png"
fig.savefig(out_plot, dpi=120)
print(f"\n[OK] {out_plot}")
