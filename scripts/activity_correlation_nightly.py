"""
RE-DO del Bloque B con datos NIGHTLY apilados (no exposiciones individuales).

El error metodológico previo: usar 1333 obs individuales infla FAP por
correlación de ruido intra-noche. Apilando a nightly:
  - N efectivo más realista
  - FAPs comparables a lo que reportamos en el resto del proyecto
  - Comparable directamente con master_v3
"""
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_CCF = "data_interim/qc/espresso_eso_ccf_curated_hr.csv"
OUTDIR    = "data_interim/qc"

BANDS = {
    'Prot 38d':  (35.0, 42.0),
    '232 d':     (225.0, 240.0),
    '270 d':     (250.0, 290.0),
    '700 d':     (650.0, 750.0),
    '53 d':      (50.0, 56.0),
    '74 d':      (60.0, 90.0),
}

# Cargar
df = pd.read_csv(INPUT_CCF)
print(f"  Cargado: N obs individuales = {len(df)}")

# Apilar a nightly por sub-segmento (model_inst)
df['bjd_int'] = df['time'].astype(int)

print("\n  Apilando a nightly (mediana por noche x segmento)...")
agg = df.groupby(['bjd_int','model_inst']).agg(
    bjd=('time','mean'),
    rv_centered=('rv_centered','median'),
    rv_err=('rv_err', lambda x: np.sqrt(np.sum(x**2))/len(x)),
    ccf_fwhm=('ccf_fwhm','median'),
    ccf_bispan=('ccf_bispan','median'),
    ccf_contrast=('ccf_contrast','median'),
    n_exp=('time','count'),
).reset_index()
print(f"  Nightly: N={len(agg)}, baseline {(agg['bjd'].max()-agg['bjd'].min())/365.25:.2f} yr")

# Centrar cada indicador (por segmento) para que la comparación sea fair
for col in ['rv_centered','ccf_fwhm','ccf_bispan','ccf_contrast']:
    agg[col + '_c'] = agg.groupby('model_inst')[col].transform(lambda x: x - x.median())

# Frecuencias
P_min, P_max = 1.5, 1500.0
nf = 30000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0/freqs

# GLS por indicador
indicators = {
    'rv':       'rv_centered_c',
    'fwhm':     'ccf_fwhm_c',
    'bispan':   'ccf_bispan_c',
    'contrast': 'ccf_contrast_c',
}

print(f"\n  GLS sobre nightly apilado...")
gls_results = {}
t = agg['bjd'].values

for label, col in indicators.items():
    y = agg[col].values
    m_ok = ~np.isnan(y)
    ls = LombScargle(t[m_ok], y[m_ok])
    pwr = ls.power(freqs)
    gls_results[label] = (pwr, ls)

# Reporte por banda
print(f"\n  {'banda':<12} {'ind':<10} {'P_pk':>8} {'pow':>7} {'FAP':>10} {'sig':>5}")
print(f"  {'-'*55}")
results = []
for bname,(plo,phi) in BANDS.items():
    in_b = (periods >= plo) & (periods <= phi)
    for label in indicators:
        pwr, ls = gls_results[label]
        i_p = np.argmax(pwr*in_b)
        P_p = periods[in_b][np.argmax(pwr[in_b])]
        pwr_p = pwr[in_b][np.argmax(pwr[in_b])]
        fap = ls.false_alarm_probability(pwr_p, method='baluev')
        sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
        print(f"  {bname:<12} {label:<10} {P_p:>7.2f} {pwr_p:>7.3f} {fap:>10.2e} {sig:>4}")
        results.append({'band':bname,'ind':label,'P_peak':P_p,'pow':pwr_p,'FAP':fap})

df_res = pd.DataFrame(results)

print(f"\n  VEREDICTO honesto por banda (CCF significativo = actividad):")
print(f"  {'banda':<12} {'rv FAP':>10} {'fwhm FAP':>10} {'bispan FAP':>10} {'contrast FAP':>12}  veredicto")
for bname in BANDS:
    sub = df_res[df_res['band']==bname]
    fap_rv = sub[sub['ind']=='rv']['FAP'].iloc[0]
    fap_f  = sub[sub['ind']=='fwhm']['FAP'].iloc[0]
    fap_b  = sub[sub['ind']=='bispan']['FAP'].iloc[0]
    fap_c  = sub[sub['ind']=='contrast']['FAP'].iloc[0]
    n_ccf_sig = sum(1 for f in [fap_f, fap_b, fap_c] if f < 0.05)
    
    if n_ccf_sig >= 2 and fap_rv < 0.05:
        verdict = f"ACTIVIDAD ({n_ccf_sig}/3 CCFs)"
    elif n_ccf_sig >= 1 and fap_rv < 0.05:
        verdict = f"actividad parcial ({n_ccf_sig}/3)"
    elif fap_rv < 0.05 and n_ccf_sig == 0:
        verdict = "RV-only → potencial Doppler"
    else:
        verdict = "sin firma significativa"
    print(f"  {bname:<12} {fap_rv:>10.2e} {fap_f:>10.2e} {fap_b:>10.2e} {fap_c:>11.2e}  {verdict}")

df_res.to_csv(f"{OUTDIR}/activity_correlation_nightly.csv", index=False)
print(f"\n  [OK] {OUTDIR}/activity_correlation_nightly.csv")

# Pearson nightly
print(f"\n  Pearson rv vs CCF (NIGHTLY):")
for ind_col in ['ccf_fwhm_c','ccf_bispan_c','ccf_contrast_c']:
    rv = agg['rv_centered_c'].values
    cc = agg[ind_col].values
    m = ~(np.isnan(rv)|np.isnan(cc))
    if m.sum()>10:
        r,p = pearsonr(rv[m],cc[m])
        sig = "***" if p<0.001 else ("**" if p<0.01 else ("*" if p<0.10 else " "))
        print(f"    rv vs {ind_col:<20} r={r:>+6.3f}  p={p:.2e} {sig}")

# Plot
fig, axes = plt.subplots(4,1,figsize=(13,10),sharex=True)
labels = {'rv':'RV','fwhm':'CCF FWHM','bispan':'CCF Bispan','contrast':'CCF Contrast'}
colors = {'rv':'black','fwhm':'tab:blue','bispan':'tab:orange','contrast':'tab:green'}

for ax,label in zip(axes, indicators.keys()):
    pwr, ls = gls_results[label]
    fap_th = ls.false_alarm_level([0.10,0.01,0.001], method='baluev')
    ax.plot(periods, pwr, color=colors[label], lw=0.7)
    for fp in fap_th:
        ax.axhline(fp, color='gray', ls='--', lw=0.5, alpha=0.5)
    for bname,(plo,phi) in BANDS.items():
        ax.axvspan(plo, phi, color='yellow', alpha=0.15)
    ax.set_xscale('log')
    ax.set_ylabel(labels[label])
    ax.grid(alpha=0.3)
axes[0].set_title('GLS rv y CCF indicators sobre NIGHTLY apilado')
axes[-1].set_xlabel('Period (d)')
fig.tight_layout()
fig.savefig(f"{OUTDIR}/activity_correlation_nightly.png", dpi=120)
print(f"  [OK] {OUTDIR}/activity_correlation_nightly.png")
