"""
Verificación final de actividad cromosférica en banda 232 d:
  - Sin detrending (centrado por mediana)
  - Detrending pol1 (lo que hicimos antes)
  - Detrending pol2

Si Hα mantiene FAP muy bajo en los tres → confirmado.
"""
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_CCF   = "data_interim/qc/espresso_eso_ccf_curated_hr.csv"
INPUT_CHROM = "data_interim/qc/activity_indicators_s1d.csv"

# Reconstruir nightly como antes
df_ccf = pd.read_csv(INPUT_CCF)
df_ccf['bjd_int'] = df_ccf['time'].astype(int)
agg = df_ccf.groupby(['bjd_int','model_inst']).agg(
    bjd=('time','mean'),
    rv=('rv_centered','median'),
).reset_index()

df_chrom = pd.read_csv(INPUT_CHROM)
df_chrom['bjd'] = df_chrom['mjd_obs'] + 2400000.5
df_chrom['bjd_int'] = df_chrom['bjd'].astype(int)
agg_chrom = df_chrom.groupby('bjd_int').agg(
    Halpha=('Halpha_idx','median'),
    CaII_K=('CaII_K_idx','median'),
).reset_index()

agg = agg.merge(agg_chrom, on='bjd_int', how='left')

t = agg['bjd'].values

P_min, P_max = 1.5, 500.0
nf = 30000
freqs = np.linspace(1.0/P_max, 1.0/P_min, nf)
periods = 1.0/freqs

BANDS = {'232 d':(225,240), '270 d':(250,290), 'Prot 38d':(35,42)}

print("="*72)
print("VERIFICACION 232/270 d en Halpha vs metodo de detrending")
print("="*72)

for col in ['rv', 'Halpha', 'CaII_K']:
    print(f"\n  Indicador: {col}")
    y_raw = agg[col].values
    m_ok = ~np.isnan(y_raw)
    y_pre = y_raw[m_ok]
    t_use = t[m_ok]
    
    # Tres versiones
    versions = {
        'sin detrend (centrar)': y_pre - np.median(y_pre),
        'detrend pol1':           y_pre - np.polyval(np.polyfit(t_use, y_pre, 1), t_use),
        'detrend pol2':           y_pre - np.polyval(np.polyfit(t_use, y_pre, 2), t_use),
    }
    
    for vname, y_v in versions.items():
        ls = LombScargle(t_use, y_v)
        pwr = ls.power(freqs)
        print(f"    [{vname}]")
        for bname,(plo,phi) in BANDS.items():
            in_b = (periods >= plo) & (periods <= phi)
            i_p = np.argmax(pwr*in_b)
            P_p = periods[in_b][np.argmax(pwr[in_b])]
            pwr_p = pwr[in_b][np.argmax(pwr[in_b])]
            fap = ls.false_alarm_probability(pwr_p, method='baluev')
            sig = "***" if fap<0.001 else ("**" if fap<0.01 else ("*" if fap<0.10 else " "))
            print(f"      {bname:<10}  P={P_p:>7.2f}  pow={pwr_p:>6.3f}  FAP={fap:.2e}  {sig}")
