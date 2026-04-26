"""
Capa 2 — Descarga RVs de Alpha Cen B desde DACE y construye masters por
sub-instrumento.

DACE devuelve dict anidado: {INST: {DRS: {MODE: {col: ndarray}}}}
Lo aplanamos a un DataFrame por sub-instrumento y guardamos como CSV.
"""

import os
import sys
import numpy as np
import pandas as pd
from dace_query.spectroscopy import Spectroscopy

OUTDIR = "data_interim/qc"
os.makedirs(OUTDIR, exist_ok=True)

# Fechas del fiber upgrade HARPS y del cambio de criostato ESPRESSO
BJD_HARPS_FIBER_START   = 2457161.5   # 16 mayo 2015
BJD_HARPS_FIBER_END     = 2457173.5   # 28 mayo 2015
BJD_ESPRESSO_CRYO_START = 2458640.0   # ~ jun 2019, intervencion ESPRESSO
BJD_ESPRESSO_CRYO_END   = 2458650.0

print("[INFO] Consultando DACE para HD128621 (Alpha Cen B)...")
data = Spectroscopy.get_timeseries(
    target='HD128621',
    sorted_by_instrument=True,
    output_format='pandas'
)

# ---------------------------------------------------------------
# Aplanar a DataFrames y construir masters por sub-instrumento
# ---------------------------------------------------------------
all_masters = []
summary = []

for inst, by_drs in data.items():
    inst_str = str(inst)
    for drs, by_mode in by_drs.items():
        drs_str = str(drs)
        for mode, col_dict in by_mode.items():
            mode_str = str(mode) if str(mode) else 'DEFAULT'

            # Reconstruir DataFrame
            try:
                df = pd.DataFrame({k: np.asarray(v) for k, v in col_dict.items()})
            except Exception as e:
                print(f"[WARN] No pude armar DF para {inst_str}/{drs_str}/{mode_str}: {e}")
                continue

            n = len(df)
            if n == 0:
                continue

            # BJD = rjd + 2400000
            df['bjd'] = df['rjd'] + 2400000.0

            # Etiqueta de sub-instrumento
            sub_inst = f"{inst_str}_{mode_str}" if mode_str != 'DEFAULT' else inst_str

            # Para HARPS metemos pre/post split de seguridad (DACE ya lo hace
            # pero lo redundamos en una columna explicita para facilitar joins)
            if 'HARPS' in inst_str:
                df['inst_split'] = np.where(
                    df['bjd'] < BJD_HARPS_FIBER_START, 'HARPS_pre',
                    np.where(df['bjd'] > BJD_HARPS_FIBER_END, 'HARPS_post', 'HARPS_GAP')
                )
            elif 'ESPRESSO' in inst_str:
                df['inst_split'] = np.where(
                    df['bjd'] < BJD_ESPRESSO_CRYO_START, 'ESPRESSO_pre',
                    np.where(df['bjd'] > BJD_ESPRESSO_CRYO_END, 'ESPRESSO_post', 'ESPRESSO_GAP')
                )
            else:
                df['inst_split'] = inst_str

            df['dace_inst'] = inst_str
            df['dace_drs']  = drs_str
            df['dace_mode'] = mode_str
            df['sub_inst']  = sub_inst

            # Stats
            bjd_min = df['bjd'].min()
            bjd_max = df['bjd'].max()
            yr = (bjd_max - bjd_min) / 365.25
            rv_med = df['rv'].median()
            rv_std = df['rv'].std()
            erv_med = df['rv_err'].median()
            qc_pass = int((df['drs_qc'] == 1).sum()) if 'drs_qc' in df else None

            summary.append({
                'sub_inst': sub_inst,
                'N': n,
                'qc_pass': qc_pass,
                'bjd_min': bjd_min,
                'bjd_max': bjd_max,
                'baseline_yr': yr,
                'rv_med': rv_med,
                'rv_std': rv_std,
                'erv_med': erv_med,
            })

            # Guardar CSV individual
            safe_name = sub_inst.replace('/', '_').replace(' ', '_')
            out = f"{OUTDIR}/dace_{safe_name}.csv"
            df.to_csv(out, index=False)
            print(f"[OK] {sub_inst}: N={n}, baseline={yr:.2f}yr, "
                  f"BJD {bjd_min:.0f}->{bjd_max:.0f} -> {out}")

            all_masters.append(df)

# ---------------------------------------------------------------
# Master unificado todos los instrumentos DACE
# ---------------------------------------------------------------
if all_masters:
    df_all = pd.concat(all_masters, ignore_index=True).sort_values('bjd').reset_index(drop=True)
    out_all = f"{OUTDIR}/dace_full_master.csv"
    df_all.to_csv(out_all, index=False)
    print(f"\n[OK] Master DACE unificado: N={len(df_all)} -> {out_all}")
    print(f"     Sub-instrumentos: {df_all['sub_inst'].value_counts().to_dict()}")

# ---------------------------------------------------------------
# Resumen tabular
# ---------------------------------------------------------------
print("\n" + "="*78)
print("RESUMEN POR SUB-INSTRUMENTO")
print("="*78)
df_sum = pd.DataFrame(summary)
if not df_sum.empty:
    df_sum = df_sum.sort_values('bjd_min').reset_index(drop=True)
    pd.set_option('display.width', 200)
    pd.set_option('display.max_columns', 20)
    print(df_sum.to_string(index=False, float_format=lambda x: f'{x:.3f}'))
    df_sum.to_csv(f"{OUTDIR}/dace_summary.csv", index=False)
    print(f"\n[OK] Resumen -> {OUTDIR}/dace_summary.csv")
