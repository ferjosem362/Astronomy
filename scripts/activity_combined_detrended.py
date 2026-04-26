"""
Análisis combinado de actividad estelar:
  1. Detrending lineal de RV y CCF nightly
  2. Apilar cromosféricos S1D (Halpha, NaD, CaII) a nightly
  3. GLS sobre 8 indicadores (4 CCF + 4 cromosféricos)
  4. Veredicto cruzado por banda

Bandas validas (>1 ciclo en baseline 1.4 yr):
  Prot 38d, 53d, 74d, 232d, 270d
"""

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from astropy.time import Time
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_CCF   = "data_interim/qc/espresso_eso_ccf_curated_hr.csv"
INPUT_CHROM = "data_interim/qc/activity_indicators_s1d.csv"
OUTDIR      = "data_interim/qc"

BANDS = {
    'Prot 38d':  (35.0, 42.0),
    '53 d':      (50.0, 56.0),
    '74 d':      (60.0, 90.0),
    '232 d':     (225.0, 240.0),
    '270 d':     (250.0, 290.0),
}

# ============================================================
# CARGAR Y APILAR CCF nightly
# ============================================================
print("="*72)
print("ANALISIS COMBINADO ACTIVIDAD: detrending + cromosfericos")
print("="*72)

df_ccf = pd.read_csv(INPUT_CCF)
print(f"\n  CCF dataset: {len(df_ccf)} obs individuales")

df_ccf['bjd_int'] = df_ccf['time'].astype(int)
agg = df_ccf.groupby(['bjd_int','model_inst']).agg(
    bjd=('time','mean'),
    rv=('rv_centered','median'),
    rv_err=('rv_err', lambda x: np.sqrt(np.sum(x**2))/len(x)),
    fwhm=('ccf_fwhm','median'),
    bispan=('ccf_bispan','median'),
    contrast=('ccf_contrast','median'),
    n_exp=('time','count'),
).reset_index()
print(f"  CCF nightly: N={len(agg)}, baseline {(agg['bjd'].max()-agg['bjd'].min())/365.25:.2f} yr")

# ============================================================
# CARGAR Y APILAR cromosfericos nightly
# ============================================================
df_chrom = pd.read_csv(INPUT_CHROM)
print(f"\n  Cromosfericos S1D: {len(df_chrom)} obs individuales")

# date_obs -> bjd (aprox via mjd_obs)
df_chrom['bjd'] = df_chrom['mjd_obs'] + 2400000.5
df_chrom['bjd_int'] = df_chrom['bjd'].astype(int)

# Apilar por noche (cromosfericos no tienen sub-segmento, tomar mediana)
chrom_cols = ['Halpha_idx', 'NaD1_idx', 'CaII_K_idx', 'CaII_H_idx']
agg_chrom = df_chrom.groupby('bjd_int').agg(
    bjd_chrom=('bjd', 'mean'),
    Halpha=('Halpha_idx', 'median'),
    NaD1=('NaD1_idx', 'median'),
    CaII_K=('CaII_K_idx', 'median'),
    CaII_H=('CaII_H_idx', 'median'),
    n_chrom=('bjd', 'count'),
).reset_index()
print(f"  Cromosfericos nightly: N={len(agg_chrom)}")

# Merge CCF nightly con cromosfericos nightly por bjd_int
agg = agg.merge(agg_chrom, on='bjd_int', how='left')
n_match = agg['Halpha'].notna().sum()
print(f"  Match CCF-Chrom: {n_match}/{len(agg)} noches con ambos")

# ============================================================
# DETRENDING lineal por sub-segmento
# ============================================================
print(f"\n  Detrending lineal por sub-segmento (model_inst)...")

indicators_to_detrend = ['rv', 'fwhm', 'bispan', 'contrast',
                         'Halpha', 'NaD1', 'CaII_K', 'CaII_H']

for col in indicators_to_detrend:
    if col not in agg.columns:
        continue
    detrended = []
    for seg in agg['model_inst'].unique():
        m = (agg['model_inst'] == seg) & agg[col].notna()
        if m.sum() < 5:
            agg.loc[m, col + '_dt'] = agg.loc[m, col] - agg.loc[m, col].mean()
            continue
        t_seg = agg.loc[m, 'bjd'].values
        y_seg = agg.loc[m, col].values
        # Polinomio orden 1 (lineal)
        coef = np.polyfit(t_seg, y_seg, 1)
        residual = y_seg - np.polyval(coef, t_seg)
        agg.loc[m, col + '_dt'] = residual

# Stats post-detrending
print(f"\n  Std post-detrending por indicador:")
for col in indicators_to_detrend:
    if col + '_dt' in agg.columns:
        v = agg[col + '_dt'].dropna()
        if len(v) > 5:
            print(f"    {col:<10} N={len(v):>4}  std_dt={v.std():.4f}  "
                  f"std_pre={agg[col].dropna().std():.4f}")

# ============================================================
# Pearson detrendeada (¿la correlación se mantiene?)
# ============================================================
print(f"\n  Pearson POST-detrending (rv vs cada indicador):")
for col in indicators_to_detrend[1:]:  # skip rv
    if col + '_dt' not in agg.columns:
        continue
    rv = agg['rv_dt'].values
    cc = agg[col + '_dt'].values
    m = ~(np.isnan(rv) | np.isnan(cc))
    if m.sum() > 10:
        r, p = pearsonr(rv[m], cc[m])
        sig = "***" if p<0.001 else ("**" if p<0.01 else ("*" if p<0.10 else " "))
        print(f"    rv_dt vs {col + '_dt':<14} r={r:>+6.3f}  p={p:.2e} {sig}  N={m.sum()}")

# ============================================================
# GLS sobre detrendeados
# ============================================================
print(f"\n  GLS sobre detrendeados...")

P_min, P_max = 1.5, 500.0  # techo 500 d (>1 ciclo en 1.4 yr)
nf = 30000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0/freqs

t = agg['bjd'].values
gls_results = {}

for col in indicators_to_detrend:
    col_dt = col + '_dt'
    if col_dt not in agg.columns:
        continue
    y = agg[col_dt].values
    m = ~np.isnan(y)
    if m.sum() < 30:
        print(f"    {col}: N={m.sum()} insuficiente, skip")
        continue
    ls = LombScargle(t[m], y[m])
    pwr = ls.power(freqs)
    gls_results[col] = (pwr, ls, m.sum())

# ============================================================
# Reporte por banda
# ============================================================
print(f"\n  REPORTE POR BANDA (detrendeado):")
print(f"  {'banda':<10} {'indicador':<10} {'N':>5} {'P_pk':>7} {'pow':>7} {'FAP':>10} {'sig':>5}")
print(f"  {'-'*60}")

results = []
for bname,(plo,phi) in BANDS.items():
    in_b = (periods >= plo) & (periods <= phi)
    for col in indicators_to_detrend:
        if col not in gls_results:
            continue
        pwr, ls, n_used = gls_results[col]
        i_p = np.argmax(pwr*in_b)
        P_p = periods[in_b][np.argmax(pwr[in_b])]
        pwr_p = pwr[in_b][np.argmax(pwr[in_b])]
        fap = ls.false_alarm_probability(pwr_p, method='baluev')
        sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
        print(f"  {bname:<10} {col:<10} {n_used:>5} {P_p:>7.2f} {pwr_p:>7.3f} {fap:>10.2e} {sig:>4}")
        results.append({'band':bname, 'ind':col, 'P_peak':P_p, 'pow':pwr_p, 'FAP':fap, 'N':n_used})

df_res = pd.DataFrame(results)

# ============================================================
# Veredicto cruzado
# ============================================================
print(f"\n  VEREDICTO CRUZADO:")
print(f"  Banda       | RV    | CCF (3) | Crom (4) | Veredicto final")
print(f"  {'-'*78}")

ccf_inds = ['fwhm', 'bispan', 'contrast']
chrom_inds = ['Halpha', 'NaD1', 'CaII_K', 'CaII_H']

for bname in BANDS:
    sub = df_res[df_res['band']==bname]
    
    rv_row = sub[sub['ind']=='rv']
    rv_fap = rv_row['FAP'].iloc[0] if len(rv_row) else 1.0
    rv_sig = "RV+" if rv_fap < 0.05 else "rv-"
    
    ccf_sub = sub[sub['ind'].isin(ccf_inds)]
    n_ccf_sig = (ccf_sub['FAP'] < 0.05).sum() if len(ccf_sub) else 0
    
    chrom_sub = sub[sub['ind'].isin(chrom_inds)]
    n_chrom_sig = (chrom_sub['FAP'] < 0.05).sum() if len(chrom_sub) else 0
    
    if rv_fap < 0.05 and (n_ccf_sig + n_chrom_sig) >= 3:
        verdict = "ACTIVIDAD CONFIRMADA (CCF+Crom)"
    elif rv_fap < 0.05 and n_chrom_sig >= 1:
        verdict = "actividad probable (Crom apoya)"
    elif rv_fap < 0.05 and n_ccf_sig >= 1:
        verdict = "actividad parcial (CCF only)"
    elif rv_fap < 0.05 and n_ccf_sig + n_chrom_sig == 0:
        verdict = "RV-only → POTENCIAL DOPPLER"
    elif rv_fap >= 0.05:
        verdict = "no firma en RV"
    else:
        verdict = "ambiguo"
    
    print(f"  {bname:<11} | {rv_sig:<5} | {n_ccf_sig}/3     | {n_chrom_sig}/4      | {verdict}")

df_res.to_csv(f"{OUTDIR}/activity_combined_detrended.csv", index=False)
print(f"\n  [OK] {OUTDIR}/activity_combined_detrended.csv")

# ============================================================
# Plot 8 paneles (4 CCF + 4 Crom)
# ============================================================
fig, axes = plt.subplots(4, 2, figsize=(14, 12), sharex=True)
axes = axes.flatten()

ind_list = ['rv', 'fwhm', 'bispan', 'contrast', 'Halpha', 'NaD1', 'CaII_K', 'CaII_H']
ind_labels = {
    'rv':'RV', 'fwhm':'CCF FWHM', 'bispan':'CCF Bispan', 'contrast':'CCF Contrast',
    'Halpha':'Hα', 'NaD1':'NaD₁', 'CaII_K':'Ca II K', 'CaII_H':'Ca II H',
}
colors_ind = {
    'rv':'black', 'fwhm':'tab:blue', 'bispan':'tab:orange', 'contrast':'tab:green',
    'Halpha':'crimson', 'NaD1':'darkviolet', 'CaII_K':'teal', 'CaII_H':'goldenrod',
}

for ax, col in zip(axes, ind_list):
    if col not in gls_results:
        ax.text(0.5, 0.5, f'{col}: insuficientes datos', ha='center', va='center',
                transform=ax.transAxes)
        continue
    pwr, ls, _ = gls_results[col]
    fap_th = ls.false_alarm_level([0.10, 0.01], method='baluev')
    ax.plot(periods, pwr, color=colors_ind[col], lw=0.7)
    for fp in fap_th:
        ax.axhline(fp, color='gray', ls='--', lw=0.5, alpha=0.5)
    for bname,(plo,phi) in BANDS.items():
        ax.axvspan(plo, phi, color='yellow', alpha=0.15)
    ax.set_xscale('log')
    ax.set_ylabel(ind_labels[col], fontsize=9)
    ax.grid(alpha=0.3)

axes[-1].set_xlabel('Period (d)')
axes[-2].set_xlabel('Period (d)')
fig.suptitle('GLS detrendeado: 4 CCF + 4 cromosfericos S1D (banda amarilla = bandas testeadas)',
             y=1.001)
fig.tight_layout()
fig.savefig(f"{OUTDIR}/activity_combined_detrended.png", dpi=120)
print(f"  [OK] {OUTDIR}/activity_combined_detrended.png")
