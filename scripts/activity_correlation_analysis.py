"""
Bloque A: Extracción de indicadores cromosféricos desde S1D ESPRESSO
  - Hα equivalent width (proxy actividad cromosférica)
  - NaD doublet (5890/5896) si está cubierto
  - Ca II H&K (3933/3968) si está cubierto

Bloque B: Correlación rv vs CCF indicators (ya tenemos en archivo A)
  Bandas: Prot ~38d, 232d, 270d, 700d (esperadas actividad)
          53d, 74d (esperadas NO actividad, validación)

Output:
  data_interim/qc/activity_indicators_s1d.csv  (cromosféricos por época)
  data_interim/qc/activity_correlation_summary.csv
  data_interim/qc/activity_correlation_diagnostic.png
"""

import os
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

WORKSPACE = Path("data_work/serval/espresso_s1d/fits")
INDEX_S1D = "data_raw/metadata/espresso_s1d_index.csv"
INPUT_CCF = "data_interim/qc/espresso_eso_ccf_curated_hr.csv"
OUTDIR    = "data_interim/qc"

# Bandas a testear
BANDS = {
    'Prot 38d':  (35.0, 42.0),
    '232 d':     (225.0, 240.0),
    '270 d':     (250.0, 290.0),
    '700 d':     (650.0, 750.0),
    '53 d':      (50.0, 56.0),
    '74 d':      (60.0, 90.0),
}

# Líneas espectrales (Å, en reposo)
LINES = {
    'Halpha': {'center': 6562.808, 'core_hw': 0.3, 'cont_blue': (6555, 6557),
               'cont_red':  (6570, 6572)},
    'NaD1':   {'center': 5895.92,  'core_hw': 0.3, 'cont_blue': (5872, 5876),
               'cont_red':  (5896, 5900)},
    'CaII_K': {'center': 3933.66,  'core_hw': 1.0, 'cont_blue': (3925, 3930),
               'cont_red':  (3937, 3942)},
    'CaII_H': {'center': 3968.47,  'core_hw': 1.0, 'cont_blue': (3960, 3965),
               'cont_red':  (3972, 3976)},
}

C_KMS = 299792.458

# ============================================================
# BLOQUE A — extraer cromosféricos de S1D
# ============================================================
print("="*72)
print("BLOQUE A — EXTRACCION CROMOSFERICA SOBRE S1D ESPRESSO")
print("="*72)

# Leer index
idx = pd.read_csv(INDEX_S1D)
print(f"\n  Index S1D: {len(idx)} entries")
# Filtrar solo HR (UHR tiene cobertura distinta)
mask_hr = idx['ins_mode'].str.contains('HR', na=False) & ~idx['ins_mode'].str.contains('UHR', na=False)
idx_hr = idx[mask_hr].copy()
print(f"  HR mode: {len(idx_hr)} entries")

# Inspeccionar 1 archivo para ver columnas reales del SPECTRUM table
sample_path = WORKSPACE / Path(idx_hr.iloc[0]['filename']).name
if not sample_path.exists():
    # buscar por nombre
    matches = list(WORKSPACE.glob(f"*{idx_hr.iloc[0]['filename']}*"))
    if matches: sample_path = matches[0]
    else:
        # último recurso: primer archivo del workspace
        all_fits = sorted(WORKSPACE.glob("*.fits"))
        sample_path = all_fits[0] if all_fits else None

if sample_path is None or not sample_path.exists():
    print(f"  [ERROR] no encuentro un FITS de muestra")
    raise SystemExit

print(f"\n  Inspecting structure: {sample_path.name}")
with fits.open(sample_path) as hdul:
    print(f"    HDU1 columns: {hdul[1].data.columns.names}")
    print(f"    HDU1 row 0 - sample data:")
    row = hdul[1].data[0]
    for col in hdul[1].data.columns.names:
        try:
            v = row[col]
            if hasattr(v, '__len__'):
                print(f"      {col}: array len={len(v)}, "
                      f"min={np.nanmin(v):.3f}, max={np.nanmax(v):.3f}")
            else:
                print(f"      {col}: {v}")
        except Exception as e:
            print(f"      {col}: error - {e}")

# Función de extracción
def extract_line_indicator(fits_path, line_info, berv_kms=0.0):
    """Extrae EW en core - continuum, corregido por BERV."""
    try:
        with fits.open(fits_path) as hdul:
            spec = hdul[1].data[0]
            wave = spec['WAVE'] if 'WAVE' in hdul[1].data.columns.names else spec['wavelength']
            flux = spec['FLUX'] if 'FLUX' in hdul[1].data.columns.names else spec['flux']
            err  = spec['ERR']  if 'ERR'  in hdul[1].data.columns.names else \
                   (spec['error'] if 'error' in hdul[1].data.columns.names else None)
        
        # Corregir BERV (datos ESO ya suelen estar en bary frame, pero por si acaso)
        # ESO S1D: WAVE ya está en sistema baricéntrico, no hace falta corregir
        wave_rest = wave  # sin corrección extra
        
        c = line_info['center']
        hw = line_info['core_hw']
        cb_lo, cb_hi = line_info['cont_blue']
        cr_lo, cr_hi = line_info['cont_red']
        
        # Verificar cobertura
        if c < wave_rest.min() or c > wave_rest.max():
            return None, None, "out_of_range"
        if cb_lo < wave_rest.min() or cr_hi > wave_rest.max():
            return None, None, "continuum_out_of_range"
        
        # Núcleo: ventana ±core_hw alrededor de la línea
        m_core = (wave_rest >= c - hw) & (wave_rest <= c + hw)
        # Continuo: media de las dos ventanas laterales
        m_cb = (wave_rest >= cb_lo) & (wave_rest <= cb_hi)
        m_cr = (wave_rest >= cr_lo) & (wave_rest <= cr_hi)
        
        if m_core.sum() < 3 or m_cb.sum() < 3 or m_cr.sum() < 3:
            return None, None, "insufficient_pixels"
        
        flux_core = np.nanmedian(flux[m_core])
        flux_cont = (np.nanmedian(flux[m_cb]) + np.nanmedian(flux[m_cr])) / 2
        
        if flux_cont <= 0:
            return None, None, "bad_continuum"
        
        index = flux_core / flux_cont  # ratio core/continuum
        # err proxy: stdev del flujo en el núcleo / N
        err_proxy = np.nanstd(flux[m_core]) / np.sqrt(m_core.sum()) / flux_cont
        
        return index, err_proxy, "ok"
    except Exception as e:
        return None, None, f"error:{type(e).__name__}"

# Procesar TODOS los HR
print(f"\n  Procesando {len(idx_hr)} archivos HR...")
records = []
status_count = {}

for i, (_, row) in enumerate(idx_hr.iterrows()):
    fname = Path(row['filename']).name if 'filename' in row else None
    if fname is None:
        continue
    fpath = WORKSPACE / fname
    if not fpath.exists():
        # buscar
        candidates = list(WORKSPACE.glob(f"*{fname}*"))
        if candidates: fpath = candidates[0]
        else:
            status_count['file_missing'] = status_count.get('file_missing', 0) + 1
            continue
    
    rec = {
        'mjd_obs': row.get('mjd_obs', np.nan),
        'date_obs': row.get('date_obs', None),
        'filename': fname,
        'ins_mode': row.get('ins_mode', None),
        'snr': row.get('snr_header', np.nan),
        'berv_kms': row.get('berv_kms', np.nan),
        'ccf_rv_kms': row.get('ccf_rv_kms', np.nan),
        'ccf_fwhm_kms': row.get('ccf_fwhm_kms', np.nan),
    }
    
    for line_name, line_info in LINES.items():
        idx_val, idx_err, status = extract_line_indicator(fpath, line_info, rec['berv_kms'])
        rec[f'{line_name}_idx'] = idx_val
        rec[f'{line_name}_err'] = idx_err
        rec[f'{line_name}_st']  = status
        status_count[f'{line_name}:{status}'] = status_count.get(f'{line_name}:{status}', 0) + 1
    
    records.append(rec)
    
    if (i+1) % 20 == 0:
        print(f"    {i+1}/{len(idx_hr)} processed")

df_chrom = pd.DataFrame(records)
out_chrom = f"{OUTDIR}/activity_indicators_s1d.csv"
df_chrom.to_csv(out_chrom, index=False)
print(f"\n  [OK] {out_chrom}: {len(df_chrom)} entries")

print(f"\n  Status por línea:")
for k, n in sorted(status_count.items()):
    print(f"    {k}: {n}")

# Para cada línea, reportar disponibilidad
print(f"\n  Resumen indicadores extraídos OK:")
for line in LINES:
    col = f'{line}_idx'
    if col in df_chrom.columns:
        n_ok = df_chrom[col].notna().sum()
        if n_ok > 10:
            v = df_chrom[col].dropna()
            print(f"    {line}: N={n_ok}, mediana={v.median():.4f}, std={v.std():.4f}")
        else:
            print(f"    {line}: N={n_ok} (insuficiente, descartar)")

# ============================================================
# BLOQUE B — correlación rv vs CCF indicators
# ============================================================
print("\n" + "="*72)
print("BLOQUE B — CORRELACION RV vs CCF INDICATORS")
print("="*72)

df_ccf = pd.read_csv(INPUT_CCF)
print(f"\n  Dataset CCF: N={len(df_ccf)}")
print(f"  Cols disponibles: {[c for c in df_ccf.columns if 'rv' in c or 'ccf' in c]}")

# Usar rv_centered (ya post-Pourbaix por mediana) y los CCF indicators
df_ccf['rv_use'] = df_ccf['rv_centered']

# Frecuencias para GLS
P_min, P_max = 1.5, 3000.0
nf = 30000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0/freqs

# Calcular GLS de cada indicador y reportar pico en cada banda
indicators = ['rv_use', 'ccf_fwhm', 'ccf_bispan', 'ccf_contrast']
gls_results = {}

t = df_ccf['time'].values
print(f"\n  Calculando GLS para cada indicador...")
for ind in indicators:
    y = df_ccf[ind].values
    if np.isnan(y).any():
        m_ok = ~np.isnan(y)
        t_i, y_i = t[m_ok], y[m_ok]
    else:
        t_i, y_i = t, y
    # detrend lineal
    p = np.polyfit(t_i, y_i, 1)
    y_dt = y_i - np.polyval(p, t_i)
    
    ls = LombScargle(t_i, y_dt)
    pwr = ls.power(freqs)
    gls_results[ind] = (t_i, y_dt, pwr, ls)

# Reporte por banda
print(f"\n  {'banda':<14} {'indicador':<14} {'P_peak':>9} {'pow':>7} {'FAP':>10} {'sig':>4}")
print(f"  {'-'*60}")
correlation_results = []
for bname, (plo, phi) in BANDS.items():
    in_b = (periods >= plo) & (periods <= phi)
    for ind in indicators:
        _, _, pwr, ls = gls_results[ind]
        i_p = np.argmax(pwr * in_b)
        P_p = periods[in_b][np.argmax(pwr[in_b])]
        pwr_p = pwr[in_b][np.argmax(pwr[in_b])]
        fap = ls.false_alarm_probability(pwr_p, method='baluev')
        sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else "-"))
        print(f"  {bname:<14} {ind:<14} {P_p:>8.2f} {pwr_p:>7.3f} {fap:>10.2e} {sig:>4}")
        correlation_results.append({
            'band': bname, 'indicator': ind, 'P_peak': P_p,
            'power': pwr_p, 'FAP': fap,
        })

df_corr = pd.DataFrame(correlation_results)
df_corr.to_csv(f"{OUTDIR}/activity_correlation_summary.csv", index=False)
print(f"\n  [OK] {OUTDIR}/activity_correlation_summary.csv")

# Veredicto banda por banda
print(f"\n  VEREDICTO POR BANDA (significativo si CCF tiene FAP<1% en banda):")
for bname in BANDS:
    sub = df_corr[df_corr['band'] == bname]
    rv_fap  = sub[sub['indicator']=='rv_use']['FAP'].iloc[0]
    ccf_faps = sub[sub['indicator']!='rv_use']['FAP'].values
    n_ccf_sig = (ccf_faps < 0.01).sum()
    
    rv_marker = "RV-sig" if rv_fap < 0.01 else "rv-no"
    if n_ccf_sig >= 2:
        verdict = f"ACTIVIDAD probable ({n_ccf_sig}/3 CCFs con FAP<1%)"
    elif n_ccf_sig == 1:
        verdict = f"actividad parcial (1/3 CCFs)"
    else:
        verdict = "sin firma de actividad → potencial Doppler"
    print(f"    {bname:<14} | {rv_marker:<8} | {verdict}")

# ============================================================
# Pearson rv vs CCF (correlación directa)
# ============================================================
print(f"\n  CORRELACION DIRECTA rv vs CCF (Pearson)")
for ind in ['ccf_fwhm','ccf_bispan','ccf_contrast']:
    rv = df_ccf['rv_use'].values
    cc = df_ccf[ind].values
    m = ~(np.isnan(rv) | np.isnan(cc))
    if m.sum() > 10:
        r, p = pearsonr(rv[m], cc[m])
        sig = "***" if p<0.001 else ("**" if p<0.01 else ("*" if p<0.10 else " "))
        print(f"    rv vs {ind:<15} r={r:>+6.3f}  p={p:.2e} {sig}")

# ============================================================
# PLOT
# ============================================================
fig, axes = plt.subplots(4, 1, figsize=(14, 12), sharex=True)

ind_labels = {'rv_use':'RV (centered)', 'ccf_fwhm':'CCF FWHM',
              'ccf_bispan':'CCF Bispan', 'ccf_contrast':'CCF Contrast'}
colors_ind = {'rv_use':'black', 'ccf_fwhm':'tab:blue',
              'ccf_bispan':'tab:orange', 'ccf_contrast':'tab:green'}

for ax, ind in zip(axes, indicators):
    _, _, pwr, ls = gls_results[ind]
    fap_th = ls.false_alarm_level([0.10, 0.01, 0.001], method='baluev')
    ax.plot(periods, pwr, color=colors_ind[ind], lw=0.7)
    for fl, fp in zip([10,1,0.1], fap_th):
        ax.axhline(fp, color='gray', ls='--', lw=0.5, alpha=0.5)
    
    for bname, (plo, phi) in BANDS.items():
        ax.axvspan(plo, phi, color='yellow', alpha=0.15)
    
    ax.set_xscale('log')
    ax.set_ylabel(ind_labels[ind])
    ax.grid(alpha=0.3)
    if ax is axes[0]:
        ax.set_title('GLS de RV y CCF indicators (banda amarilla = bandas testeadas)')

axes[-1].set_xlabel('Period (days)')
fig.tight_layout()
out = f"{OUTDIR}/activity_correlation_diagnostic.png"
fig.savefig(out, dpi=120)
print(f"\n  [OK] {out}")

print(f"\n{'='*72}")
print("DONE — interpretación con datos en mano")
print('='*72)
