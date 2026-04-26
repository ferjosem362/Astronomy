"""
Build master_v3: master_v2 + ESPRESSO 2023-2024 expandido (1333 obs).

FASE 1: Aplicar Pourbaix2016 al archivo curated_hr (RV crudas).
FASE 2: Filtrar outliers, apilar a nightly.
FASE 3: Mergear a master_v2 -> master_v3 (sin pisar ESPRESSO19_HR11).
FASE 4: GLS combinado sobre master_v3, comparar con v2.
"""

import os
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_EXT     = "data_interim/qc/espresso_eso_ccf_curated_hr.csv"
INPUT_MASTER  = "data_interim/qc/master_v2_nightly_residuals_clean.csv"
OUTDIR        = "data_interim/qc"

OUT_NIGHTLY   = f"{OUTDIR}/espresso_external_nightly.csv"
OUT_MASTER_V3 = f"{OUTDIR}/master_v3_nightly_residuals_clean.csv"
OUT_PLOT      = f"{OUTDIR}/master_v3_diagnostic.png"

# Pourbaix2016 (mismos parametros que master_v2)
P_yr_orb = 79.91
e_orb    = 0.524
omega_B  = np.deg2rad(232.3)
T0_BJD   = 2435349.815
V0_ms    = -22390.0
dVB      = 329.0
K_B      = 5503.97  # m/s

# ============================================================
# FASE 1 — Pourbaix sobre archivo externo
# ============================================================
print(f"\n{'='*72}")
print("FASE 1: POURBAIX2016 SOBRE ESPRESSO EXTERNO")
print('='*72)

df = pd.read_csv(INPUT_EXT)
print(f"\n  N={len(df)}, BJD {df['time'].min():.0f}-{df['time'].max():.0f} "
      f"({(df['time'].max()-df['time'].min())/365.25:.2f} yr)")
print(f"  Sub-segmentos (model_inst):")
for seg in sorted(df['model_inst'].unique()):
    print(f"    {seg}: N={(df['model_inst']==seg).sum()}")

def true_anomaly(t, P_d, e, T0):
    M = 2*np.pi*(t - T0)/P_d % (2*np.pi)
    E = M.copy()
    for _ in range(150):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < 1e-13: break
    nu = 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    return nu

nu = true_anomaly(df['time'].values, P_yr_orb*365.25, e_orb, T0_BJD)
rv_model = V0_ms + K_B * (e_orb*np.cos(omega_B) + np.cos(omega_B + nu)) + dVB
df['rv_resid'] = df['rv'] - rv_model
df['rv_resid_centered'] = df.groupby('model_inst')['rv_resid'].transform(
    lambda x: x - x.median())

print(f"\n  rv mediana cruda:     {df['rv'].median():.2f} m/s")
print(f"  rv_model mediana:     {pd.Series(rv_model).median():.2f} m/s")
print(f"\n  Std post-Pourbaix por segmento (centrado):")
for seg in sorted(df['model_inst'].unique()):
    m = df['model_inst']==seg
    sub = df.loc[m, 'rv_resid_centered']
    print(f"    {seg:<28} N={m.sum():>4}  std={sub.std():.2f} m/s")

# ============================================================
# FASE 2 — Filtro outliers + nightly
# ============================================================
print(f"\n{'='*72}")
print("FASE 2: FILTRO OUTLIERS + NIGHTLY APILADO")
print('='*72)

def robust_sigma(x):
    return 1.4826 * (x - x.median()).abs().median()

df['sigma_r'] = df.groupby('model_inst')['rv_resid_centered'].transform(robust_sigma)
df['is_outlier'] = df['rv_resid_centered'].abs() > 5 * df['sigma_r']
n_out = df['is_outlier'].sum()
print(f"\n  Outliers MAD 5sigma: {n_out}/{len(df)} ({n_out/len(df)*100:.1f}%)")

df_clean = df[~df['is_outlier']].copy()
print(f"  N tras filtro: {len(df_clean)}")
for seg in sorted(df_clean['model_inst'].unique()):
    m = df_clean['model_inst']==seg
    print(f"    {seg:<28} N={m.sum():>4}  "
          f"std_clean={df_clean.loc[m,'rv_resid_centered'].std():.3f} m/s")

# Apilar a nightly (mediana por noche por segmento)
df_clean['bjd_int'] = df_clean['time'].astype(int)
agg = df_clean.groupby(['bjd_int','model_inst']).agg(
    bjd=('time','mean'),
    rv=('rv','mean'),
    rv_resid=('rv_resid_centered','mean'),
    rv_err=('rv_err', lambda x: np.sqrt(np.sum(x**2))/len(x)),
    n_exp=('time','count'),
    ccf_fwhm=('ccf_fwhm','mean'),
    ccf_bispan=('ccf_bispan','mean'),
    ccf_contrast=('ccf_contrast','mean'),
).reset_index()
agg = agg.rename(columns={'model_inst':'sub_inst_orig'})

# Renombrar para no chocar con ESPRESSO19_HR11 del master_v2
seg_map = {
    'ESPRESSO_HR_2023_EARLY': 'ESPRESSO_EXT_2023e',
    'ESPRESSO_HR_2023PLUS':   'ESPRESSO_EXT_2023plus',
}
agg['sub_inst'] = agg['sub_inst_orig'].map(seg_map).fillna(agg['sub_inst_orig'])
agg = agg.sort_values('bjd').reset_index(drop=True)

print(f"\n  Nightly apilado: {len(agg)} noches")
for seg in sorted(agg['sub_inst'].unique()):
    m = agg['sub_inst']==seg
    print(f"    {seg:<28} N_nights={m.sum():>4}  "
          f"std_resid={agg.loc[m,'rv_resid'].std():.3f} m/s")

agg.to_csv(OUT_NIGHTLY, index=False)
print(f"\n  [OK] {OUT_NIGHTLY}")

# ============================================================
# FASE 3 — Merge a master_v3
# ============================================================
print(f"\n{'='*72}")
print("FASE 3: MERGE A MASTER_V3")
print('='*72)

master_v2 = pd.read_csv(INPUT_MASTER)
print(f"\n  Master_v2: {len(master_v2)} noches")
print(f"  Sub-instrumentos master_v2:")
for sub in sorted(master_v2['sub_inst'].unique()):
    n = (master_v2['sub_inst']==sub).sum()
    print(f"    {sub:<28} N={n}")

# Cross-check con ESPRESSO19_HR11 actual: ¿hay solapamiento?
esp19 = master_v2[master_v2['sub_inst']=='ESPRESSO19_HR11'].copy()
esp19['bjd_int'] = esp19['bjd'].astype(int)
agg['bjd_int'] = agg['bjd'].astype(int)

ovr = pd.merge(
    esp19[['bjd_int','rv_resid','rv_err']].rename(
        columns={'rv_resid':'rv_v2','rv_err':'err_v2'}),
    agg[['bjd_int','rv_resid','rv_err','sub_inst']].rename(
        columns={'rv_resid':'rv_ext','rv_err':'err_ext'}),
    on='bjd_int', how='inner')

print(f"\n  Solapamiento ESPRESSO19_HR11(v2) vs externo: {len(ovr)} noches")
if len(ovr) > 0:
    diff = ovr['rv_v2'] - ovr['rv_ext']
    print(f"    Diferencia v2-ext:  mediana={diff.median():+.3f} m/s  "
          f"std={diff.std():.3f} m/s  rango=[{diff.min():+.2f},{diff.max():+.2f}]")
    if abs(diff.median()) > 2.0:
        print(f"    ⚠ OFFSET >2 m/s — anyado igual pero anotado")
    else:
        print(f"    ✓ Pipelines consistentes (offset < 2 m/s)")

# Conservar master_v2 ESPRESSO19_HR11; añadir solo nights del externo no solapadas
ovr_bjds = set(ovr['bjd_int'])
agg_to_add = agg[~agg['bjd_int'].isin(ovr_bjds)].copy()
print(f"\n  Noches del externo a anyadir (no solapadas): {len(agg_to_add)}")

# Alinear columnas para concat
master_cols = list(master_v2.columns)
ext_cols = ['bjd','sub_inst','rv','rv_err','rv_resid']

new_rows = pd.DataFrame()
for col in master_cols:
    if col in agg_to_add.columns:
        new_rows[col] = agg_to_add[col].values
    else:
        new_rows[col] = np.nan

# Añadir CCF indicators como columnas extras (master_v2 tendrá NaN)
for col in ['ccf_fwhm','ccf_bispan','ccf_contrast','n_exp']:
    if col not in master_v2.columns:
        master_v2[col] = np.nan
    new_rows[col] = agg_to_add[col].values

master_v3 = pd.concat([master_v2, new_rows], ignore_index=True)
master_v3 = master_v3.sort_values('bjd').reset_index(drop=True)

print(f"\n  Master_v3: {len(master_v3)} noches "
      f"(v2={len(master_v2)} + {len(new_rows)} ESPRESSO 2023-2024)")
print(f"\n  Sub-instrumentos master_v3:")
for sub in sorted(master_v3['sub_inst'].unique()):
    m = master_v3['sub_inst']==sub
    std = master_v3.loc[m,'rv_resid'].std()
    print(f"    {sub:<28} N={m.sum():>4}  std_resid={std:.3f}")

master_v3.to_csv(OUT_MASTER_V3, index=False)
print(f"\n  [OK] {OUT_MASTER_V3}")

# ============================================================
# FASE 4 — GLS combinado v2 vs v3
# ============================================================
print(f"\n{'='*72}")
print("FASE 4: GLS COMBINADO v2 vs v3 (comparativo)")
print('='*72)

# v3
master_v3['rv_resid_centered'] = master_v3.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
t3 = master_v3['bjd'].values
y3 = master_v3['rv_resid_centered'].values
dy3 = np.maximum(master_v3['rv_err'].values, 0.1)

# v2
master_v2_re = pd.read_csv(INPUT_MASTER)
master_v2_re['rv_resid_centered'] = master_v2_re.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())
t2 = master_v2_re['bjd'].values
y2 = master_v2_re['rv_resid_centered'].values
dy2 = np.maximum(master_v2_re['rv_err'].values, 0.1)

P_min, P_max = 1.5, 5000.0
nf = 50000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0/freqs

ls2 = LombScargle(t2, y2, dy2, normalization='standard')
power_v2 = ls2.power(freqs)
ls3 = LombScargle(t3, y3, dy3, normalization='standard')
power_v3 = ls3.power(freqs)

# Bandas (incluyendo ~454 d nueva, del pipeline paralelo)
TARGET_BANDS = {
    'Prot ~38 d': (35.0, 42.0),
    '~74 d':      (60.0, 90.0),
    '~234 d':     (225.0, 245.0),  # banda externa
    '~270 d':     (250.0, 290.0),
    '~454 d':     (440.0, 470.0),  # banda externa
    '~700 d':     (650.0, 750.0),
}

print(f"\n  {'banda':<14} | {'v2 P':>8} {'v2 FAP':>10} | {'v3 P':>8} {'v3 FAP':>10} | {'cambio':>10}")
print(f"  {'-'*14} | {'-'*8} {'-'*10} | {'-'*8} {'-'*10} | {'-'*10}")
for bname,(plo,phi) in TARGET_BANDS.items():
    in_b = (periods >= plo) & (periods <= phi)
    if in_b.sum() == 0: continue
    
    iv2 = np.argmax(power_v2*in_b); P2 = periods[iv2]; pw2 = power_v2[iv2]
    fap2 = ls2.false_alarm_probability(pw2, method='baluev')
    
    iv3 = np.argmax(power_v3*in_b); P3 = periods[iv3]; pw3 = power_v3[iv3]
    fap3 = ls3.false_alarm_probability(pw3, method='baluev')
    
    if fap3 < 0.01 and fap2 >= 0.01: cambio = "EMERGE"
    elif fap3 >= 0.01 and fap2 < 0.01: cambio = "DESAPARECE"
    elif fap3 < fap2/10: cambio = "REFUERZA"
    elif fap3 > fap2*10: cambio = "DEBILITA"
    else: cambio = "estable"
    
    print(f"  {bname:<14} | {P2:>7.2f} {fap2:>10.2e} | {P3:>7.2f} {fap3:>10.2e} | {cambio:>10}")

# Top 10 picos v3
peaks_idx, _ = find_peaks(power_v3, distance=50, prominence=0.005)
peaks_sorted = sorted(zip(peaks_idx, power_v3[peaks_idx]), key=lambda x:-x[1])[:10]
print(f"\n  Top 10 picos globales v3:")
for i,(idx,p) in enumerate(peaks_sorted):
    P_pk = periods[idx]
    fap = ls3.false_alarm_probability(p, method='baluev')
    sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
    print(f"    #{i+1:2d} P={P_pk:>9.3f}  pow={p:.4f}  FAP={fap:.2e} {sig}")

# Plot
fig, axes = plt.subplots(2,1,figsize=(15,9))

axes[0].plot(periods, power_v2, color='gray', lw=0.6, alpha=0.7, label=f'v2 (N={len(master_v2_re)})')
axes[0].plot(periods, power_v3, color='red', lw=0.7, alpha=0.85, label=f'v3 (N={len(master_v3)})')
fap_th2 = ls2.false_alarm_level([0.10,0.01,0.001], method='baluev')
fap_th3 = ls3.false_alarm_level([0.10,0.01,0.001], method='baluev')
for fp in fap_th3:
    axes[0].axhline(fp, color='red', ls='--', lw=0.4, alpha=0.4)
ymax = max(power_v2.max(), power_v3.max())
for bname,(plo,phi) in TARGET_BANDS.items():
    axes[0].axvspan(plo, phi, color='yellow', alpha=0.15)
    axes[0].text(np.sqrt(plo*phi), ymax*0.95, bname, ha='center', fontsize=7, alpha=0.6)
axes[0].set_xscale('log')
axes[0].set_ylabel('GLS power')
axes[0].set_title('Comparacion GLS combinado: master_v2 (gris) vs master_v3 (rojo)')
axes[0].legend(loc='upper right', fontsize=9)
axes[0].grid(alpha=0.3)

# Cobertura temporal
COLORS = {
    'HARPS03':'tab:blue','HARPS15':'tab:cyan',
    'ES_zhao2017':'tab:green','CHIRON_zhao2017':'tab:olive',
    'ESPRESSO18_HR11':'tab:red','ESPRESSO18_HR21':'tab:orange',
    'ESPRESSO18_UHR':'tab:pink','ESPRESSO19_HR11':'tab:purple',
    'ESPRESSO_EXT_2023e':'magenta','ESPRESSO_EXT_2023plus':'darkviolet',
}
for sub in sorted(master_v3['sub_inst'].unique()):
    m = master_v3['sub_inst']==sub
    c = COLORS.get(sub, 'gray')
    axes[1].errorbar(master_v3.loc[m,'bjd'], master_v3.loc[m,'rv_resid_centered'],
                     yerr=master_v3.loc[m,'rv_err'], fmt='.', ms=3, alpha=0.6,
                     color=c, label=f'{sub} N={m.sum()}')
axes[1].axhline(0, color='black', lw=0.5)
axes[1].set_xlabel('BJD')
axes[1].set_ylabel('rv_resid centrado (m/s)')
axes[1].set_title('Cobertura temporal master_v3')
axes[1].set_ylim(-15,15)
axes[1].legend(fontsize=7, loc='upper right', ncol=2)
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUT_PLOT, dpi=120)
print(f"\n  [OK] {OUT_PLOT}")

print(f"\n{'='*72}")
print("RESUMEN")
print('='*72)
print(f"  Master_v2 -> Master_v3: {len(master_v2_re)} -> {len(master_v3)} noches "
      f"(+{len(new_rows)})")
print(f"  Output: {OUT_MASTER_V3}")
print(f"  Plot:   {OUT_PLOT}")
print(f"\n  Bandas de interes que cambian → mira tabla 'cambio' arriba")
print('='*72)
