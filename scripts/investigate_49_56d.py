"""
Investigacion forense de los picos a ~49 d y ~56 d.

Tres tareas:
  A) Localizacion precisa de los picos en zoom
  B) sBGLS especifico a esas P para ver coherencia
  C) Test de aliases:
     - vs Prot=38.65 d (multiplos racionales)
     - vs window function (aliases de sampling)
     - vs sub-instrumentos individuales (consistencia)

Output:
  data_interim/qc/picos_49_56_summary.csv
  data_interim/qc/picos_49_56_diagnostic.png
"""

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"

# Banda de busqueda inicial
P_BAND_LO, P_BAND_HI = 45.0, 60.0

# Periodos candidatos a investigar (refinar despues)
P_CANDIDATES_INITIAL = [49.0, 56.0]

# Prot Alpha Cen B (de los GLS previos)
P_ROT = 38.65

# Frecuencias para zoom alta resolucion
P_min_zoom, P_max_zoom = 40.0, 65.0
nf_zoom = 30000
freqs_zoom = np.linspace(1.0/P_max_zoom, 1.0/P_min_zoom, nf_zoom)
periods_zoom = 1.0 / freqs_zoom

# Frecuencias completas (para FAP global)
nf_full = 50000
freqs_full = np.linspace(1.0/5000.0, 1.0/1.5, nf_full)
periods_full = 1.0 / freqs_full

# Cargar
df = pd.read_csv(INPUT)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
df = df.sort_values('bjd').reset_index(drop=True)

t = df['bjd'].values
y = df['rv_resid_centered'].values
dy = np.maximum(df['rv_err'].values, 0.1)

print(f"[INFO] N={len(df)}, baseline={(t.max()-t.min())/365.25:.2f} yr")

# ============================================================
# A) Localizacion precisa
# ============================================================
print("\n" + "="*70)
print("A) LOCALIZACION PRECISA EN BANDA 45-60 d")
print("="*70)

ls = LombScargle(t, y, dy, normalization='standard')
power_zoom = ls.power(freqs_zoom)

# find_peaks en zoom
peak_idx, props = find_peaks(power_zoom, distance=200, prominence=0.005)
peaks_in_band = []
for idx in peak_idx:
    P_peak = periods_zoom[idx]
    if P_BAND_LO <= P_peak <= P_BAND_HI:
        pwr = power_zoom[idx]
        fap = ls.false_alarm_probability(pwr, method='baluev')
        peaks_in_band.append((P_peak, pwr, fap))

# Ordenar por power
peaks_in_band.sort(key=lambda x: -x[1])
print(f"\n  Picos detectados en banda {P_BAND_LO}-{P_BAND_HI} d (ordenados por power):")
for i, (P, pwr, fap) in enumerate(peaks_in_band[:5]):
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"    #{i+1}  P={P:.4f} d   power={pwr:.4f}   FAP={fap:.2e}  {sig}")

# Tomar los dos picos mas altos como candidatos refinados
if len(peaks_in_band) >= 2:
    P_CAND_1 = peaks_in_band[0][0]
    P_CAND_2 = peaks_in_band[1][0]
elif len(peaks_in_band) == 1:
    P_CAND_1 = peaks_in_band[0][0]
    P_CAND_2 = None
else:
    print("  [WARN] Sin picos detectados, usando candidatos iniciales")
    P_CAND_1, P_CAND_2 = P_CANDIDATES_INITIAL

print(f"\n  Candidatos refinados: P1={P_CAND_1:.3f} d, "
      f"P2={'(none)' if P_CAND_2 is None else f'{P_CAND_2:.3f} d'}")

# ============================================================
# B) sBGLS para cada candidato
# ============================================================
print("\n" + "="*70)
print("B) sBGLS PARA CADA CANDIDATO")
print("="*70)

# Bandas estrechas para cada candidato (+/- 2%)
BANDS_REFINED = {}
for label, P in [('P1', P_CAND_1), ('P2', P_CAND_2)]:
    if P is None:
        continue
    BANDS_REFINED[f'{label}_{P:.2f}d'] = (P*0.98, P*1.02)

# Apilar cronologicamente
N_steps = list(range(50, len(df)+1, 25))
sbgls_results = {label: {'N': [], 'P_peak': [], 'power': [], 'logBF': []}
                 for label in BANDS_REFINED}

for n in N_steps:
    sub = df.iloc[:n]
    t_n = sub['bjd'].values
    y_n = sub['rv_resid_centered'].values
    dy_n = np.maximum(sub['rv_err'].values, 0.1)
    if (t_n.max() - t_n.min()) < 100:
        continue
    
    ls_n = LombScargle(t_n, y_n, dy_n, normalization='standard')
    power_n = ls_n.power(freqs_full)
    
    for label, (p_lo, p_hi) in BANDS_REFINED.items():
        in_b = (periods_full >= p_lo) & (periods_full <= p_hi)
        if in_b.sum() == 0:
            continue
        i_peak = np.argmax(power_n * in_b)
        sbgls_results[label]['N'].append(n)
        sbgls_results[label]['P_peak'].append(periods_full[i_peak])
        sbgls_results[label]['power'].append(power_n[i_peak])
        sbgls_results[label]['logBF'].append(n/2 * power_n[i_peak])

# Reporte sBGLS
for label, data in sbgls_results.items():
    if not data['N']:
        continue
    Ns = np.array(data['N'])
    Ps = np.array(data['P_peak'])
    logBF = np.array(data['logBF'])
    slope_lnbf = np.polyfit(Ns, logBF, 1)[0]
    p_relstd = Ps.std() / Ps.mean() * 100
    
    if slope_lnbf > 0.005 and p_relstd < 3:
        verdict = "COHERENTE: sigue Kepleriano"
    elif slope_lnbf > 0.005 and p_relstd > 5:
        verdict = "SOSPECHOSO: ln(BF) crece pero P migra"
    elif slope_lnbf < -0.005:
        verdict = "DECRECIENTE: actividad transitoria"
    else:
        verdict = "PLANO: sin coherencia"
    
    print(f"\n  [{label}]")
    print(f"    ln(BF) {logBF[0]:.2f} -> {logBF[-1]:.2f}, slope={slope_lnbf:.4f}")
    print(f"    P_peak: media={Ps.mean():.3f} d, std={Ps.std():.3f} d, "
          f"relstd={p_relstd:.1f}%")
    print(f"    Veredicto: {verdict}")

# ============================================================
# C) Tests de aliases
# ============================================================
print("\n" + "="*70)
print("C) TESTS DE ALIASES")
print("="*70)

# C1) vs Prot=38.65 d
print(f"\n  C1) Relacion con P_rot={P_ROT} d:")
for label, P in [('P1', P_CAND_1), ('P2', P_CAND_2)]:
    if P is None:
        continue
    ratio = P / P_ROT
    # Buscar ratio racional simple
    for n in range(1, 5):
        for m in range(1, 5):
            r_test = n/m
            if abs(ratio - r_test) < 0.05:
                print(f"    {label}={P:.3f} d / P_rot={P_ROT} = {ratio:.3f} ≈ {n}/{m}={r_test:.3f}  "
                      f"-> POSIBLE ARMONICO")
                break

# C2) vs window function
print(f"\n  C2) Window function en bandas candidatos:")
ls_win = LombScargle(t, np.ones_like(t), 1.0,
                     normalization='standard', center_data=False, fit_mean=False)
power_win = ls_win.power(freqs_full)
for label, P in [('P1', P_CAND_1), ('P2', P_CAND_2)]:
    if P is None:
        continue
    in_b = (periods_full >= P*0.95) & (periods_full <= P*1.05)
    pwr_win_max = power_win[in_b].max() if in_b.any() else 0
    print(f"    {label}={P:.3f} d: window power max = {pwr_win_max:.4f}  "
          f"({'ALIAS LIKELY' if pwr_win_max > 0.05 else 'window clean'})")

# C3) Por sub-instrumento
print(f"\n  C3) Pico por sub-instrumento (bandas estrechas):")
SUBS_BIG = ['HARPS03', 'ES_zhao2017', 'CHIRON_zhao2017', 'HARPS15', 'ESPRESSO19_HR11']
sub_results = {label: {} for label in BANDS_REFINED}

for sub_name in SUBS_BIG:
    m = df['sub_inst'] == sub_name
    if m.sum() < 10:
        continue
    t_s = df.loc[m, 'bjd'].values
    y_s = df.loc[m, 'rv_resid_centered'].values
    dy_s = np.maximum(df.loc[m, 'rv_err'].values, 0.1)
    
    ls_s = LombScargle(t_s, y_s, dy_s, normalization='standard')
    pwr_s = ls_s.power(freqs_full)
    
    for label, (p_lo, p_hi) in BANDS_REFINED.items():
        in_b = (periods_full >= p_lo) & (periods_full <= p_hi)
        if in_b.sum() == 0:
            continue
        i_peak = np.argmax(pwr_s * in_b)
        P_peak = periods_full[i_peak]
        pwr_peak = pwr_s[i_peak]
        fap_peak = ls_s.false_alarm_probability(pwr_peak, method='baluev')
        sub_results[label][sub_name] = (P_peak, pwr_peak, fap_peak)

print(f"\n    {'sub_inst':<20} ", end='')
for label in BANDS_REFINED:
    print(f"{label:>20}", end='  ')
print()

for sub_name in SUBS_BIG:
    print(f"    {sub_name:<20} ", end='')
    for label in BANDS_REFINED:
        if sub_name in sub_results[label]:
            P_p, pwr_p, fap_p = sub_results[label][sub_name]
            sig = "***" if fap_p<0.001 else ("**" if fap_p<0.01 else ("*" if fap_p<0.10 else "-"))
            print(f"  P={P_p:5.2f} fap={fap_p:.0e} {sig:>3}", end='  ')
        else:
            print(f"  {'(N<10)':>20}", end='  ')
    print()

# ============================================================
# Plot
# ============================================================
fig, axes = plt.subplots(3, 2, figsize=(15, 14))

# Top row: zoom espectrales
for ax, (label, P) in zip(axes[0], [('P1', P_CAND_1), ('P2', P_CAND_2)]):
    if P is None:
        ax.axis('off')
        continue
    P_lo_view = max(P_BAND_LO, P-5)
    P_hi_view = min(P_BAND_HI, P+5)
    mask = (periods_zoom >= P_lo_view) & (periods_zoom <= P_hi_view)
    ax.plot(periods_zoom[mask], power_zoom[mask], 'k-', lw=1)
    fap_th = ls.false_alarm_level([0.10, 0.01, 0.001], method='baluev')
    ax.axhline(fap_th[0], color='gray', ls='--', lw=0.5, alpha=0.5, label='FAP 10%')
    ax.axhline(fap_th[1], color='orange', ls='--', lw=0.6, label='FAP 1%')
    ax.axhline(fap_th[2], color='red', ls='--', lw=0.6, label='FAP 0.1%')
    ax.axvline(P, color='green', ls=':', lw=2, label=f'P={P:.2f} d')
    ax.set_xlabel('Period (d)')
    ax.set_ylabel('GLS power')
    ax.set_title(f'Zoom espectral: {label}={P:.3f} d')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

# Middle row: sBGLS
for ax, (label, data) in zip(axes[1], sbgls_results.items()):
    if not data['N']:
        ax.axis('off')
        continue
    Ns = np.array(data['N'])
    logBF = np.array(data['logBF'])
    Ps = np.array(data['P_peak'])
    
    color1 = 'tab:blue'
    ax.plot(Ns, logBF, '-', color=color1, lw=2, label='ln(BF) proxy')
    ax.fill_between(Ns, 0, logBF, alpha=0.15, color=color1)
    ax.set_ylabel('ln(BF) proxy', color=color1)
    
    ax2 = ax.twinx()
    ax2.plot(Ns, Ps, 'r.', ms=5, alpha=0.7)
    ax2.set_ylabel('P_peak (d)', color='red')
    
    ax.set_title(f'sBGLS: {label}')
    ax.set_xlabel('N epocas')
    ax.grid(alpha=0.3)

# Bottom row: comparacion sub-instrumentos
for ax, (label, P) in zip(axes[2], [('P1', P_CAND_1), ('P2', P_CAND_2)]):
    if P is None or label not in BANDS_REFINED:
        ax.axis('off')
        continue
    
    label_full = next(k for k in BANDS_REFINED if k.startswith(label))
    sub_data = sub_results.get(label_full, {})
    
    sub_names = []
    sub_powers = []
    sub_faps = []
    for sub_name in SUBS_BIG:
        if sub_name in sub_data:
            P_p, pwr_p, fap_p = sub_data[sub_name]
            sub_names.append(sub_name.replace('_zhao2017', '_z').replace('_HR11', ''))
            sub_powers.append(pwr_p)
            sub_faps.append(fap_p)
    
    if sub_names:
        colors_bar = ['green' if f < 0.01 else ('orange' if f < 0.10 else 'red')
                      for f in sub_faps]
        ax.bar(sub_names, sub_powers, color=colors_bar, edgecolor='black')
        ax.set_ylabel('GLS power en banda')
        ax.set_title(f'{label}={P:.3f} d - power por sub-inst (verde=FAP<1%)')
        ax.tick_params(axis='x', rotation=30)
        ax.grid(alpha=0.3, axis='y')

fig.tight_layout()
out = f"{OUTDIR}/picos_49_56_diagnostic.png"
fig.savefig(out, dpi=120)
print(f"\n[OK] {out}")
