"""
Stacked Bayesian Generalized Lomb-Scargle (Mortier & Cameron 2017).

Para cada banda de interes, calcula como evoluciona la potencia GLS
y la "evidencia" (proxy) en funcion del numero de epocas usadas (apiladas
cronologicamente, no en bloques aleatorios).

- Senales reales (Keplerianos): power y P_peak ESTABLES, evidencia crece monotonicamente
- Actividad transitoria: P_peak migra, evidencia sube y baja
- Aliases: P_peak salta, no estable

Salida:
  data_interim/qc/sbgls_per_band.png
  data_interim/qc/sbgls_summary.csv
"""

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"

TARGET_BANDS = {
    'Prot ~38 d':  (35.0, 42.0),
    '~74 d':       (60.0, 90.0),
    '~270 d':      (250.0, 290.0),
    '~700 d':      (650.0, 750.0),
}

df = pd.read_csv(INPUT)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
df = df.sort_values('bjd').reset_index(drop=True)

print(f"[INFO] Total epocas: {len(df)}")
print(f"[INFO] Rango temporal: {(df['bjd'].max()-df['bjd'].min())/365.25:.2f} yr")

# Apilamos cronologicamente: N=50, 75, 100, ..., total
N_steps = list(range(50, len(df)+1, 25))
if N_steps[-1] != len(df):
    N_steps.append(len(df))
print(f"[INFO] Pasos de N: {len(N_steps)} (de {N_steps[0]} a {N_steps[-1]})")

# Frecuencias
P_min, P_max_total = 1.5, 5000.0
nf = 30000
freqs = np.linspace(1.0/P_max_total, 1.0/P_min, nf)
periods = 1.0 / freqs

results = {band: {'N': [], 'P_peak': [], 'power_peak': [],
                  'logBF_proxy': [], 'baseline_yr': [], 'fap': []}
           for band in TARGET_BANDS}

for n in N_steps:
    sub_df = df.iloc[:n]
    t = sub_df['bjd'].values
    y = sub_df['rv_resid_centered'].values
    dy = np.maximum(sub_df['rv_err'].values, 0.1)
    baseline_yr = (t.max() - t.min()) / 365.25
    
    if len(t) < 10 or baseline_yr < 0.1:
        continue
    
    ls = LombScargle(t, y, dy, normalization='standard')
    power = ls.power(freqs)
    
    for band, (p_lo, p_hi) in TARGET_BANDS.items():
        in_band = (periods >= p_lo) & (periods <= p_hi)
        if in_band.sum() == 0:
            continue
        i_peak = np.argmax(power * in_band)
        p_peak = periods[i_peak]
        pow_peak = power[i_peak]
        # FAP del pico
        try:
            fap = ls.false_alarm_probability(pow_peak, method='baluev')
        except Exception:
            fap = np.nan
        # Proxy ln(BF) ~ N/2 * power (Mortier & Cameron approximation)
        logBF = n / 2.0 * pow_peak
        results[band]['N'].append(n)
        results[band]['P_peak'].append(p_peak)
        results[band]['power_peak'].append(pow_peak)
        results[band]['logBF_proxy'].append(logBF)
        results[band]['baseline_yr'].append(baseline_yr)
        results[band]['fap'].append(fap)

# ============================================================
# Plot: 4 paneles, uno por banda, 3 trazas (P_peak, power, logBF)
# ============================================================
fig, axes = plt.subplots(4, 1, figsize=(13, 14), sharex=True)

for ax, (band, data) in zip(axes, results.items()):
    if not data['N']:
        ax.set_title(f"{band}: sin datos suficientes")
        continue
    Ns = np.array(data['N'])
    pwr = np.array(data['power_peak'])
    logBF = np.array(data['logBF_proxy'])
    Ps = np.array(data['P_peak'])
    
    # Eje principal: ln(BF) proxy
    color1 = 'tab:blue'
    ax.plot(Ns, logBF, '-', color=color1, lw=2, label='ln(BF) proxy = N*power/2')
    ax.fill_between(Ns, 0, logBF, alpha=0.15, color=color1)
    ax.set_ylabel('ln(BF) proxy', color=color1)
    ax.tick_params(axis='y', labelcolor=color1)
    
    # Eje secundario: P_peak
    ax2 = ax.twinx()
    color2 = 'tab:red'
    p_lo, p_hi = TARGET_BANDS[band]
    ax2.plot(Ns, Ps, 'o', color=color2, ms=4, alpha=0.7, label='P_peak')
    ax2.axhline(p_lo, color='gray', ls=':', lw=0.5)
    ax2.axhline(p_hi, color='gray', ls=':', lw=0.5)
    ax2.set_ylabel('P_peak (d)', color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.set_ylim(p_lo - (p_hi-p_lo)*0.3, p_hi + (p_hi-p_lo)*0.3)
    
    ax.set_title(f'sBGLS: {band}')
    ax.grid(alpha=0.3)
    
    # Interpretacion
    if len(logBF) > 5:
        # Tendencia de ln(BF)
        slope_lnbf = np.polyfit(Ns, logBF, 1)[0]
        # Estabilidad de P_peak (std relativo)
        p_relstd = Ps.std() / Ps.mean() * 100
        
        if slope_lnbf > 0.005 and p_relstd < 3:
            verdict = "✓ COHERENTE: ln(BF) crece, P_peak estable"
            color_box = 'lightgreen'
        elif slope_lnbf > 0.005 and p_relstd > 5:
            verdict = "⚠ ln(BF) crece pero P_peak migra (sospechoso)"
            color_box = 'lightyellow'
        elif slope_lnbf < -0.005:
            verdict = "✗ ln(BF) decrece (transitoria/artefacto)"
            color_box = 'lightcoral'
        else:
            verdict = "~ ln(BF) plano"
            color_box = 'lightgray'
        
        ax.text(0.02, 0.95,
                f"{verdict}\nslope_lnBF={slope_lnbf:.4f}\nP_peak std={p_relstd:.1f}%",
                transform=ax.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor=color_box, alpha=0.85))

axes[-1].set_xlabel('N epocas (apiladas cronologicamente)')
fig.suptitle('Stacked Bayesian GLS: evolucion de la evidencia con N creciente', y=1.001)
fig.tight_layout()
out = f"{OUTDIR}/sbgls_per_band.png"
fig.savefig(out, dpi=120)
print(f"\n[OK] {out}")

# Salvar tabla
all_rows = []
for band, data in results.items():
    for i in range(len(data['N'])):
        all_rows.append({
            'band': band,
            'N': data['N'][i],
            'baseline_yr': data['baseline_yr'][i],
            'P_peak': data['P_peak'][i],
            'power_peak': data['power_peak'][i],
            'logBF_proxy': data['logBF_proxy'][i],
            'fap': data['fap'][i],
        })
df_summary = pd.DataFrame(all_rows)
df_summary.to_csv(f"{OUTDIR}/sbgls_summary.csv", index=False)
print(f"[OK] {OUTDIR}/sbgls_summary.csv")

# Reporte texto
print("\n" + "="*70)
print("INTERPRETACION POR BANDA")
print("="*70)
for band, data in results.items():
    if not data['N']:
        continue
    Ns = np.array(data['N'])
    Ps = np.array(data['P_peak'])
    logBF = np.array(data['logBF_proxy'])
    print(f"\n{band}:")
    print(f"  N: {Ns[0]} -> {Ns[-1]}")
    print(f"  Baseline: {data['baseline_yr'][0]:.2f} -> {data['baseline_yr'][-1]:.2f} yr")
    print(f"  P_peak: media={Ps.mean():.2f} d, std={Ps.std():.2f} d, "
          f"rango=[{Ps.min():.1f}, {Ps.max():.1f}]")
    print(f"  ln(BF): inicial={logBF[0]:.2f}, final={logBF[-1]:.2f}, "
          f"slope={np.polyfit(Ns, logBF, 1)[0]:.4f}")
