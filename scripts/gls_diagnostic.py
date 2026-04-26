"""
Fase 3 paso 2 - Periodograma GLS sobre residuos post-binario LIMPIOS.

Genera:
  - GLS por sub-instrumento (cada uno individual)
  - GLS combinado (todos los residuos centrados juntos)
  - Window function (para identificar aliases)
  - Tabla de picos significativos por encima de FAP 10%, 1%, 0.1%

Salida:
  data_interim/qc/gls_per_subinst.png
  data_interim/qc/gls_combined.png
  data_interim/qc/gls_peaks.csv
"""

import os
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"

# Bandas de interes (periodos en dias)
TARGET_BANDS = {
    '~10.9 d':    (10.0, 11.5),
    'Prot ~38 d': (35.0, 42.0),
    '~74 d':      (60.0, 90.0),
    '~270 d':     (250.0, 290.0),
    '~700 d':     (650.0, 750.0),
}

# Frecuencias a explorar: 1.5 d a 5000 d
P_min, P_max = 1.5, 5000.0
nf = 50000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0 / freqs

# Niveles FAP a marcar
FAP_LEVELS = [0.10, 0.01, 0.001]

# Colores
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

# Cargar
df = pd.read_csv(INPUT)
print(f"[INFO] Cargado {INPUT}: {len(df)} epocas")

# Centrar residuos por sub-instrumento (segunda vez por si)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())

# ============================================================
# GLS por sub-instrumento (los que tienen >=10 puntos)
# ============================================================
print("\n" + "="*70)
print("GLS POR SUB-INSTRUMENTO")
print("="*70)

valid_subs = [s for s in ORDER if (df['sub_inst']==s).sum() >= 10]
fig, axes = plt.subplots(len(valid_subs), 1, figsize=(13, 2.2*len(valid_subs)),
                         sharex=True)
if len(valid_subs) == 1:
    axes = [axes]
ax_idx = 0
peaks_table = []

for sub in ORDER:
    m = df['sub_inst'] == sub
    if m.sum() < 10:
        print(f"  [SKIP] {sub}: solo {m.sum()} puntos")
        continue
    
    t = df.loc[m, 'bjd'].values
    y = df.loc[m, 'rv_resid_centered'].values
    dy = df.loc[m, 'rv_err'].values
    dy = np.maximum(dy, 0.1)
    
    ls = LombScargle(t, y, dy, normalization='standard', center_data=True)
    power = ls.power(freqs)
    fap_thresh = ls.false_alarm_level(FAP_LEVELS, method='baluev')
    
    print(f"\n[{sub}] N={m.sum()}, baseline={(t.max()-t.min())/365.25:.2f} yr")
    
    # Picos por banda
    for band_name, (p_lo, p_hi) in TARGET_BANDS.items():
        in_band = (periods >= p_lo) & (periods <= p_hi)
        if in_band.sum() == 0:
            continue
        i_peak = np.argmax(power * in_band)
        p_peak = periods[i_peak]
        pow_peak = power[i_peak]
        fap = ls.false_alarm_probability(pow_peak, method='baluev')
        sig = "***" if fap < 0.001 else ("**" if fap < 0.01 else ("*" if fap < 0.10 else " "))
        print(f"    {band_name:15s}  P={p_peak:8.3f} d  pow={pow_peak:.4f}  FAP={fap:.2e}  {sig}")
        peaks_table.append({
            'sub_inst': sub,
            'band': band_name,
            'P_peak': p_peak,
            'power': pow_peak,
            'FAP': fap,
            'N': m.sum(),
            'baseline_yr': (t.max()-t.min())/365.25,
        })
    
    ax = axes[ax_idx]
    ax.plot(periods, power, color=COLORS[sub], lw=0.7, alpha=0.9)
    for fl, fp in zip(FAP_LEVELS, fap_thresh):
        ax.axhline(fp, color='gray', ls='--', lw=0.5, alpha=0.6)
    for band_name, (p_lo, p_hi) in TARGET_BANDS.items():
        ax.axvspan(p_lo, p_hi, color='yellow', alpha=0.15)
    ax.set_xscale('log')
    ax.set_ylabel(f'{sub}\n(N={m.sum()})', fontsize=8)
    ax.grid(alpha=0.3)
    ax_idx += 1

axes[-1].set_xlabel('Period (days)')
fig.suptitle('GLS por sub-instrumento (residuos centrados, post-Pourbaix2016)', y=1.001)
fig.tight_layout()
out = f"{OUTDIR}/gls_per_subinst.png"
fig.savefig(out, dpi=120)
print(f"\n[OK] {out}")

# ============================================================
# GLS COMBINADO (todos los puntos centrados juntos)
# ============================================================
print("\n" + "="*70)
print("GLS COMBINADO (todos los sub-instrumentos)")
print("="*70)

t_all = df['bjd'].values
y_all = df['rv_resid_centered'].values
dy_all = np.maximum(df['rv_err'].values, 0.1)

ls_all = LombScargle(t_all, y_all, dy_all, normalization='standard', center_data=True)
power_all = ls_all.power(freqs)
fap_all = ls_all.false_alarm_level(FAP_LEVELS, method='baluev')

# Window function (sin restar el dato, solo el sampling)
ls_win = LombScargle(t_all, np.ones_like(t_all), 1.0,
                     normalization='standard', center_data=False, fit_mean=False)
power_win = ls_win.power(freqs)

fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)

axes[0].plot(periods, power_all, 'k-', lw=0.7)
for fl, fp in zip(FAP_LEVELS, fap_all):
    axes[0].axhline(fp, color='gray', ls='--', lw=0.6, alpha=0.7,
                    label=f'FAP {fl*100:.1f}%')
ymax = power_all.max()
for band_name, (p_lo, p_hi) in TARGET_BANDS.items():
    axes[0].axvspan(p_lo, p_hi, color='yellow', alpha=0.2)
    axes[0].text(np.sqrt(p_lo*p_hi), ymax*0.95,
                 band_name, ha='center', fontsize=8, alpha=0.7)
axes[0].set_xscale('log')
axes[0].set_ylabel('GLS power')
axes[0].set_title(f'GLS combinado (N={len(t_all)} epocas, todos los sub-inst, residuos centrados)')
axes[0].legend(loc='upper right', fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].plot(periods, power_win, 'darkred', lw=0.7)
axes[1].axvline(365.25, color='blue', ls=':', lw=1, alpha=0.7, label='1 year')
axes[1].axvline(29.53, color='green', ls=':', lw=1, alpha=0.7, label='Lunar month')
axes[1].axvline(1.0, color='orange', ls=':', lw=1, alpha=0.7, label='1 day')
for band_name, (p_lo, p_hi) in TARGET_BANDS.items():
    axes[1].axvspan(p_lo, p_hi, color='yellow', alpha=0.2)
axes[1].set_xscale('log')
axes[1].set_xlabel('Period (days)')
axes[1].set_ylabel('Window power')
axes[1].set_title('Window function (sampling) - picos aqui = aliases potenciales')
axes[1].legend(loc='upper right', fontsize=8)
axes[1].grid(alpha=0.3)

fig.tight_layout()
out2 = f"{OUTDIR}/gls_combined.png"
fig.savefig(out2, dpi=120)
print(f"[OK] {out2}")

# ============================================================
# Picos en bandas para el combinado
# ============================================================
print("\n" + "="*70)
print("PICOS EN BANDAS DE INTERES (combinado)")
print("="*70)
for band_name, (p_lo, p_hi) in TARGET_BANDS.items():
    in_band = (periods >= p_lo) & (periods <= p_hi)
    i_peak = np.argmax(power_all * in_band)
    p_peak = periods[i_peak]
    pow_peak = power_all[i_peak]
    fap = ls_all.false_alarm_probability(pow_peak, method='baluev')
    sig = "***" if fap < 0.001 else ("**" if fap < 0.01 else ("*" if fap < 0.10 else " "))
    print(f"  {band_name:15s}  P={p_peak:9.3f} d   pow={pow_peak:.4f}   FAP={fap:.2e}   {sig}")
    peaks_table.append({
        'sub_inst': 'COMBINED',
        'band': band_name,
        'P_peak': p_peak,
        'power': pow_peak,
        'FAP': fap,
        'N': len(t_all),
        'baseline_yr': (t_all.max()-t_all.min())/365.25,
    })

# Top 10 picos globales
print("\n" + "="*70)
print("TOP 10 PICOS GLOBALES (combinado)")
print("="*70)
peaks_idx, props = find_peaks(power_all, distance=50, prominence=0.005)
peaks_sorted = sorted(zip(peaks_idx, power_all[peaks_idx]),
                      key=lambda x: -x[1])[:10]
for i, (idx, p) in enumerate(peaks_sorted):
    p_peak = periods[idx]
    fap = ls_all.false_alarm_probability(p, method='baluev')
    sig = "***" if fap < 0.001 else ("**" if fap < 0.01 else ("*" if fap < 0.10 else " "))
    print(f"  #{i+1:2d}  P={p_peak:9.3f} d   pow={p:.4f}   FAP={fap:.2e}   {sig}")

df_peaks = pd.DataFrame(peaks_table)
df_peaks.to_csv(f"{OUTDIR}/gls_peaks.csv", index=False)
print(f"\n[OK] {OUTDIR}/gls_peaks.csv")

print("\n" + "="*70)
print("Leyenda: * = FAP<10%, ** = FAP<1%, *** = FAP<0.1%")
print("Bandera roja: si pico a 270d aparece tambien en window function -> alias")
print("="*70)
