"""
Cross-validacion bundle SERVAL vs master_v3_v2 en noches solapadas.

Mide consistencia entre dos extracciones distintas de las mismas
observaciones ESPRESSO 2023-2024.
"""
import sys, os
sys.path.insert(0, '.')

import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Pourbaix2016 - copiado del script original
P_yr = 79.91
P_s = P_yr * 365.25 * 86400.0
P_days = P_yr * 365.25
e_orb = 0.524
i_rad = np.deg2rad(79.32)
omega_B = np.deg2rad(232.3)
T_jyear = 1955.66
V0_ms = -22390.0
dVB_ms = 329.0

T0_BJD = Time(T_jyear, format='jyear', scale='tdb').jd
G = 6.67430e-11
M_sun = 1.98892e30
M_A = 1.133
M_B = 0.972
K_B = (2*np.pi*G * (M_A * M_sun)**3 / (P_s * ((M_A + M_B)*M_sun)**2))**(1/3) \
      * np.sin(i_rad) / np.sqrt(1 - e_orb**2)

def true_anomaly(t):
    M = 2*np.pi * (t - T0_BJD) / P_days
    M = np.mod(M, 2*np.pi)
    E = M.copy()
    for _ in range(100):
        dE = (E - e_orb*np.sin(E) - M) / (1 - e_orb*np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < 1e-12: break
    return 2 * np.arctan2(np.sqrt(1+e_orb)*np.sin(E/2), np.sqrt(1-e_orb)*np.cos(E/2))

def pourbaix_rv(bjd):
    nu = true_anomaly(bjd)
    return V0_ms + K_B * (e_orb*np.cos(omega_B) + np.cos(omega_B + nu)) + dVB_ms

# Cargar
bundle = pd.read_csv("release/espresso_sflA_independent_SERVAL_bundle/data/espresso_sflA_drs_nightly_group_detrended_clean.csv")
master = pd.read_csv("data_interim/qc/master_v3_v2_nightly_residuals_clean.csv")

# Filtrar master a ESPRESSO 2023-2024 (los dos sub_inst que cubren el bundle)
m_e = master[master.sub_inst.isin(['ESPRESSO19_HR11', 'ESPRESSO_EXT_2023plus'])].copy()
m_e_in_range = m_e[(m_e.bjd >= bundle.bjd.min()-1) & (m_e.bjd <= bundle.bjd.max()+1)].copy()

print("="*72)
print("CROSS-VALIDATION BUNDLE vs MASTER")
print("="*72)
print(f"\nBundle:                     N={len(bundle)}")
print(f"Master ESPRESSO 2023-2024:  N={len(m_e_in_range)}")
print(f"  - ESPRESSO19_HR11:        N={(m_e_in_range.sub_inst=='ESPRESSO19_HR11').sum()}")
print(f"  - ESPRESSO_EXT_2023plus:  N={(m_e_in_range.sub_inst=='ESPRESSO_EXT_2023plus').sum()}")

# Calcular Pourbaix del bundle
bundle['rv_pourbaix'] = pourbaix_rv(bundle.bjd.values)
bundle['rv_resid_bundle'] = bundle.rv_rel_ms - bundle.rv_pourbaix

# Match each bundle night with closest master night within 0.5d
matches = []
for _, r in bundle.iterrows():
    diffs = np.abs(m_e_in_range.bjd.values - r.bjd)
    if diffs.min() < 0.5:
        idx = diffs.argmin()
        match_row = m_e_in_range.iloc[idx]
        matches.append({
            'bjd': r.bjd,
            'sub_inst_master': match_row.sub_inst,
            'rv_resid_bundle_raw': r.rv_resid_bundle,
            'rv_resid_master': match_row.rv_resid,
            'rv_err_bundle': r.erv_ms,
            'sn55_bundle': r.sn55_median,
        })

if not matches:
    print("\n[ERROR] No matches found")
    sys.exit(1)

df = pd.DataFrame(matches)
print(f"\nMatched nights:  N={len(df)}")
print(f"  vs ESPRESSO19_HR11:        N={(df.sub_inst_master=='ESPRESSO19_HR11').sum()}")
print(f"  vs ESPRESSO_EXT_2023plus:  N={(df.sub_inst_master=='ESPRESSO_EXT_2023plus').sum()}")

# Calcular zero-point offset y diferencias
zero_point = df.rv_resid_bundle_raw.median() - df.rv_resid_master.median()
df['rv_resid_bundle_aligned'] = df.rv_resid_bundle_raw - zero_point
df['diff_ms'] = df.rv_resid_bundle_aligned - df.rv_resid_master

print(f"\nZero-point offset: {zero_point:+.3f} m/s (esto es V0 sistémico, esperado)")

# Stats
print(f"\n{'='*72}")
print(f"DIFERENCIAS bundle - master (despues de zero-point):")
print(f"{'='*72}")
print(f"  Mediana:  {df.diff_ms.median():+.3f} m/s")
print(f"  Media:    {df.diff_ms.mean():+.3f} m/s")
print(f"  Std:      {df.diff_ms.std():.3f} m/s")
print(f"  Min/Max:  {df.diff_ms.min():+.2f} / {df.diff_ms.max():+.2f} m/s")
print(f"  P16/P84:  {df.diff_ms.quantile(0.16):+.2f} / {df.diff_ms.quantile(0.84):+.2f}")

# Analisis por sub-bloque
for sub in ['ESPRESSO19_HR11', 'ESPRESSO_EXT_2023plus']:
    mask = df.sub_inst_master == sub
    if mask.sum() > 0:
        sub_diff = df[mask].diff_ms
        print(f"\n  {sub} (N={mask.sum()}):")
        print(f"    Mediana: {sub_diff.median():+.3f}")
        print(f"    Std:     {sub_diff.std():.3f}")

# Test trend con BJD
from scipy.stats import linregress
lr = linregress(df.bjd, df.diff_ms)
print(f"\n  Trend con BJD:")
print(f"    slope:   {lr.slope*365.25*1000:+.3f} mm/s/yr (esperado ~0)")
print(f"    p-value: {lr.pvalue:.4f}")
if lr.pvalue < 0.01:
    print(f"    [WARN] Trend significativo - extracciones drift entre s\u00ed")
else:
    print(f"    [OK] Sin trend significativo")

# Test trend con SNR
lr_sn = linregress(df.sn55_bundle, df.diff_ms)
print(f"\n  Trend con SNR:")
print(f"    slope:   {lr_sn.slope:+.4f} m/s per SNR unit")
print(f"    p-value: {lr_sn.pvalue:.4f}")
if lr_sn.pvalue < 0.01:
    print(f"    [WARN] Diferencia depende de SNR")
else:
    print(f"    [OK] Sin dependencia con SNR")

# Save
OUT = "data_interim/qc/bundle_master_crossval.csv"
df.to_csv(OUT, index=False)
print(f"\n[OK] {OUT}")

# Plot
fig, axes = plt.subplots(2, 2, figsize=(14, 9))

# Scatter master vs bundle
ax = axes[0,0]
ax.errorbar(df.rv_resid_master, df.rv_resid_bundle_aligned,
            xerr=2.0, yerr=df.rv_err_bundle, fmt='o', alpha=0.6, ms=5)
xmin, xmax = df.rv_resid_master.min(), df.rv_resid_master.max()
ax.plot([xmin, xmax], [xmin, xmax], 'k--', alpha=0.5, label='1:1')
ax.set_xlabel('rv_resid master (m/s)')
ax.set_ylabel('rv_resid bundle aligned (m/s)')
ax.set_title(f'Master vs Bundle: N={len(df)}, std diff={df.diff_ms.std():.2f} m/s')
ax.legend()
ax.grid(alpha=0.3)

# Diff vs BJD
ax = axes[0,1]
for sub, color in [('ESPRESSO19_HR11', 'tab:purple'), ('ESPRESSO_EXT_2023plus', 'tab:orange')]:
    mask = df.sub_inst_master == sub
    if mask.sum() > 0:
        ax.scatter(df[mask].bjd, df[mask].diff_ms, color=color, s=20, alpha=0.7,
                   label=f'{sub} (N={mask.sum()})')
ax.axhline(0, color='k', lw=0.5)
ax.axhline(df.diff_ms.median(), color='red', ls='--', alpha=0.5,
           label=f'median={df.diff_ms.median():+.2f}')
ax.set_xlabel('BJD')
ax.set_ylabel('bundle - master (m/s)')
ax.set_title(f'Trend vs BJD (slope p={lr.pvalue:.3f})')
ax.legend(fontsize=8)
ax.grid(alpha=0.3)

# Histograma diff
ax = axes[1,0]
ax.hist(df.diff_ms, bins=20, alpha=0.7, edgecolor='k')
ax.axvline(0, color='k', lw=0.5)
ax.axvline(df.diff_ms.median(), color='red', ls='--',
           label=f'median={df.diff_ms.median():+.2f}')
ax.set_xlabel('bundle - master (m/s)')
ax.set_ylabel('N')
ax.set_title(f'Histograma diferencias (std={df.diff_ms.std():.2f})')
ax.legend()
ax.grid(alpha=0.3)

# Diff vs SNR
ax = axes[1,1]
ax.scatter(df.sn55_bundle, df.diff_ms, s=20, alpha=0.7)
ax.axhline(0, color='k', lw=0.5)
ax.set_xlabel('SNR55 bundle')
ax.set_ylabel('bundle - master (m/s)')
ax.set_title(f'Trend vs SNR (p={lr_sn.pvalue:.3f})')
ax.grid(alpha=0.3)

fig.suptitle('Cross-validation: Bundle SERVAL vs Master DACE')
fig.tight_layout()
fig.savefig("data_interim/qc/bundle_master_crossval.png", dpi=120)
print(f"[OK] data_interim/qc/bundle_master_crossval.png")

# Veredicto final
print(f"\n{'='*72}")
print(f"VEREDICTO")
print(f"{'='*72}")
std = df.diff_ms.std()
if std < 1.0:
    print(f"  [EXCELENTE] std diff = {std:.2f} m/s < 1 m/s")
    print(f"  Las dos extracciones (DACE master vs SERVAL bundle) son CONSISTENTES.")
    print(f"  Paper puede citar cross-validacion como sanity check fuerte.")
elif std < 3.0:
    print(f"  [BUENO] std diff = {std:.2f} m/s entre 1-3 m/s")
    print(f"  Diferencias compatibles con jitter ESPRESSO normal entre extracciones.")
    print(f"  Paper documenta como diferencia metodologica menor.")
else:
    print(f"  [ATENCION] std diff = {std:.2f} m/s > 3 m/s")
    print(f"  Hay diferencia metodologica seria entre las dos extracciones.")
    print(f"  Paper debe explicar el origen de la discrepancia.")
