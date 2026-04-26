"""
Fase 3 paso 3 - Inject-recovery test.

Inyecta planetas Keplerianos circulares a:
  - P=53 d con K ∈ {0.0317 (1.19xMars), 0.1, 0.3, 0.5, 1.0, 2.0, 5.0} m/s
  - P=74 d con K ∈ {0.0388 (1.65xMars), 0.1, 0.3, 0.5, 1.0, 2.0, 5.0} m/s

Para cada (P, K):
  1. Suma senial Kepleriana en epocas reales del master nightly
  2. Recalcula GLS combinado
  3. Reporta: power en banda P_inj +/- 5%, FAP, recuperabilidad

Output:
  data_interim/qc/inject_recovery_results.csv
  data_interim/qc/inject_recovery_diagnostic.png
"""

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR = "data_interim/qc"

# Parametros estelares
M_STAR_SUN = 0.972
M_MARS_EARTH = 0.107   # M_Marte en Tierras

# Periodos a testear
PERIODS = [53.0, 74.0]

# Amplitudes K en m/s. Incluye 'tu' caso (Mars-equivalente) y grilla de mapeo
K_GRID = [
    0.0317,  # 1.19 x Mars en M_planet → K(P=53)
    0.0388,  # 1.65 x Mars en M_planet → K(P=74)
    0.1, 0.3, 0.5, 1.0, 2.0, 5.0
]

# Bandas para reportar
BAND_HALF_WIDTH_REL = 0.05  # ±5% del periodo inyectado

# Frecuencias GLS
P_min, P_max = 1.5, 5000.0
nf = 50000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods_grid = 1.0 / freqs

# Cargar
df = pd.read_csv(INPUT)
df['rv_resid_centered'] = df.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
print(f"[INFO] Master clean: N={len(df)}, baseline={(df['bjd'].max()-df['bjd'].min())/365.25:.2f} yr")

# Orden temporal
df = df.sort_values('bjd').reset_index(drop=True)

t = df['bjd'].values
y_orig = df['rv_resid_centered'].values
dy = np.maximum(df['rv_err'].values, 0.1)

# Calcular FAP del baseline (sin inyeccion)
ls_base = LombScargle(t, y_orig, dy, normalization='standard')
power_base = ls_base.power(freqs)

# Picos en bandas de los inyectados (BASELINE, sin inyeccion)
print("\n" + "="*70)
print("BASELINE (sin inyeccion) - power en bandas de interes")
print("="*70)
baseline_band_power = {}
for P_inj in PERIODS:
    p_lo = P_inj * (1 - BAND_HALF_WIDTH_REL)
    p_hi = P_inj * (1 + BAND_HALF_WIDTH_REL)
    in_band = (periods_grid >= p_lo) & (periods_grid <= p_hi)
    i_peak = np.argmax(power_base * in_band)
    P_peak = periods_grid[i_peak]
    pow_peak = power_base[i_peak]
    fap = ls_base.false_alarm_probability(pow_peak, method='baluev')
    baseline_band_power[P_inj] = pow_peak
    print(f"  P_inj={P_inj} d → P_peak_band={P_peak:.3f} d, power={pow_peak:.4f}, FAP={fap:.2e}")

# ============================================================
# Inyeccion + recuperacion
# ============================================================
results = []

print("\n" + "="*70)
print("INYECCION-RECUPERACION")
print("="*70)
print(f"{'P_inj':>6} {'K_inj':>8} {'M·sini':>10} "
      f"{'P_peak':>9} {'pow_peak':>9} {'pow_excess':>11} "
      f"{'FAP':>10} {'recovered':>10}")
print("-"*78)

for P_inj in PERIODS:
    p_lo = P_inj * (1 - BAND_HALF_WIDTH_REL)
    p_hi = P_inj * (1 + BAND_HALF_WIDTH_REL)
    in_band = (periods_grid >= p_lo) & (periods_grid <= p_hi)
    
    for K_inj in K_GRID:
        # Senial Kepleriana circular: K * sin(2π t/P + phi0)
        # phi0 = 0 (no importa para deteccion)
        signal = K_inj * np.sin(2 * np.pi * t / P_inj)
        y_inj = y_orig + signal
        
        # GLS sobre y_inj
        ls = LombScargle(t, y_inj, dy, normalization='standard')
        power = ls.power(freqs)
        
        # Pico en banda inyectada
        i_peak = np.argmax(power * in_band)
        P_peak = periods_grid[i_peak]
        pow_peak = power[i_peak]
        pow_excess = pow_peak - baseline_band_power[P_inj]
        fap = ls.false_alarm_probability(pow_peak, method='baluev')
        
        # M_planeta en M_Tierra desde K
        # K = 0.6395 * (M_p sini / M_⊕) * (M_*/M_sun)^(-2/3) * (P/yr)^(-1/3)  [m/s]
        # Despejando: M_p sini = K / 0.6395 * (M_*)^(2/3) * (P/yr)^(1/3)
        P_yr = P_inj / 365.25
        M_p_earth = K_inj / 0.6395 * (M_STAR_SUN)**(2/3) * (P_yr)**(1/3)
        M_p_mars = M_p_earth / M_MARS_EARTH
        
        # Criterio de recuperacion: power_inj > power_baseline + 0.02 (5sigma sobre baseline)
        # Y FAP < 1%
        recovered = (pow_excess > 0.02) and (fap < 0.01) and (abs(P_peak - P_inj) / P_inj < 0.03)
        rec_str = "YES" if recovered else "no"
        
        print(f"{P_inj:>6.0f} {K_inj:>8.4f} {M_p_earth:>7.3f}M⊕ "
              f"{P_peak:>9.3f} {pow_peak:>9.4f} {pow_excess:>+11.4f} "
              f"{fap:>10.2e} {rec_str:>10}")
        
        results.append({
            'P_inj': P_inj,
            'K_inj': K_inj,
            'M_p_earth': M_p_earth,
            'M_p_mars': M_p_mars,
            'P_peak': P_peak,
            'power_peak': pow_peak,
            'power_baseline': baseline_band_power[P_inj],
            'power_excess': pow_excess,
            'FAP': fap,
            'recovered': recovered,
        })
    print()

df_res = pd.DataFrame(results)
df_res.to_csv(f"{OUTDIR}/inject_recovery_results.csv", index=False)
print(f"[OK] {OUTDIR}/inject_recovery_results.csv")

# ============================================================
# Plot recuperabilidad vs K
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

for ax, P_inj in zip(axes, PERIODS):
    sub = df_res[df_res['P_inj'] == P_inj]
    
    rec = sub[sub['recovered']]
    not_rec = sub[~sub['recovered']]
    
    ax.scatter(rec['K_inj'], rec['power_excess'], s=80, color='tab:green',
               label='Recovered', zorder=3, edgecolor='black')
    ax.scatter(not_rec['K_inj'], not_rec['power_excess'], s=80, color='tab:red',
               label='Not recovered', zorder=3, edgecolor='black', alpha=0.7)
    
    # Anotar masas
    for _, r in sub.iterrows():
        ax.annotate(f"{r['M_p_earth']:.2f}M⊕", 
                    xy=(r['K_inj'], r['power_excess']),
                    xytext=(5, 5), textcoords='offset points', fontsize=7)
    
    ax.set_xscale('log')
    ax.set_yscale('symlog', linthresh=0.001)
    ax.axhline(0.02, color='gray', ls='--', lw=1, label='Detection threshold (excess=0.02)')
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xlabel('K_inj (m/s)')
    if ax is axes[0]:
        ax.set_ylabel('Power excess over baseline')
    ax.set_title(f'P_inj = {P_inj} d')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3, which='both')

fig.suptitle('Inject-recovery: amplitud detectable a 53 d y 74 d')
fig.tight_layout()
out = f"{OUTDIR}/inject_recovery_diagnostic.png"
fig.savefig(out, dpi=120)
print(f"[OK] {out}")

# ============================================================
# Resumen: K minima detectable
# ============================================================
print("\n" + "="*70)
print("LIMITES DE DETECCION")
print("="*70)
for P_inj in PERIODS:
    sub = df_res[df_res['P_inj'] == P_inj].sort_values('K_inj')
    rec = sub[sub['recovered']]
    if len(rec) > 0:
        K_min = rec['K_inj'].min()
        M_min = rec.loc[rec['K_inj'].idxmin(), 'M_p_earth']
        print(f"  P={P_inj} d: K_min recuperada = {K_min:.3f} m/s "
              f"({M_min:.2f} M⊕, {M_min/M_MARS_EARTH:.1f} x Marte)")
    else:
        print(f"  P={P_inj} d: ninguna K recuperada en grilla {K_GRID}")
        print(f"    -> upper limit pendiente de extender la grilla a K mayor")

print("="*70)
