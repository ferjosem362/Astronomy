"""
QC diagnostico de los CSV DACE recien bajados.

Investiga:
1. is_latest_drs en HARPS03 (porque hay 19343 puntos cuando esperabamos ~2000)
2. distribucion temporal HARPS03 (anos antes de 2008?)
3. outliers en ESPRESSO19 (std=1857 es enorme)
4. drs_qc por sub-instrumento
"""

import pandas as pd
import numpy as np

QC_FILES = {
    'HARPS03':       'data_interim/qc/dace_HARPS03_HARPS.csv',
    'HARPS15':       'data_interim/qc/dace_HARPS15_HARPS.csv',
    'ESPRESSO18_HR11': 'data_interim/qc/dace_ESPRESSO18_SINGLEHR11.csv',
    'ESPRESSO18_HR21': 'data_interim/qc/dace_ESPRESSO18_SINGLEHR21.csv',
    'ESPRESSO18_UHR':  'data_interim/qc/dace_ESPRESSO18_SINGLEUHR.csv',
    'ESPRESSO19_HR11': 'data_interim/qc/dace_ESPRESSO19_SINGLEHR11.csv',
}

for name, path in QC_FILES.items():
    print("="*78)
    print(f"  {name}")
    print("="*78)
    df = pd.read_csv(path)
    print(f"N total: {len(df)}")
    
    # is_latest_drs
    if 'is_latest_drs' in df.columns:
        print(f"\nis_latest_drs:")
        print(df['is_latest_drs'].value_counts())
    
    # drs_qc
    if 'drs_qc' in df.columns:
        print(f"\ndrs_qc:")
        print(df['drs_qc'].value_counts())
    
    # Distribucion temporal por anio
    df['year'] = pd.to_datetime(df['bjd'] - 2400000.5, unit='D', origin='1858-11-17').dt.year
    print(f"\nN por anio:")
    print(df['year'].value_counts().sort_index())
    
    # Stats con QC=1 + is_latest_drs=True
    if 'is_latest_drs' in df.columns and 'drs_qc' in df.columns:
        clean = df[(df['is_latest_drs'] == True) & (df['drs_qc'] == 1)]
        print(f"\nDespues de filtros (is_latest_drs=True AND drs_qc=1): N={len(clean)}")
        if len(clean) > 0:
            print(f"  RV mediana: {clean['rv'].median():.4f}")
            print(f"  RV std:     {clean['rv'].std():.4f}")
            print(f"  eRV mediana: {clean['rv_err'].median():.4f}")
            print(f"  RV min:  {clean['rv'].min():.2f}")
            print(f"  RV max:  {clean['rv'].max():.2f}")
            
            # detectar outliers > 5 sigma
            rv_med = clean['rv'].median()
            rv_mad = (clean['rv'] - rv_med).abs().median()
            outlier_thresh = 5 * 1.4826 * rv_mad
            n_out = ((clean['rv'] - rv_med).abs() > outlier_thresh).sum()
            print(f"  Outliers > 5 MAD ({outlier_thresh:.1f} m/s del med): {n_out}")
    
    # Programas / pub_reference
    if 'pub_reference' in df.columns:
        refs = df['pub_reference'].dropna().value_counts().head(5)
        if len(refs) > 0:
            print(f"\nTop pub_references:")
            print(refs)
    
    if 'program_id' in df.columns:
        progs = df['program_id'].dropna().value_counts().head(5)
        if len(progs) > 0:
            print(f"\nTop program_id:")
            print(progs)
    
    print()
