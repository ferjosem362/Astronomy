"""
Diagnostico rapido de los CSV principales UVES y CHIRON existentes,
para decidir como adoptarlos al esquema DACE.
"""
import os
import pandas as pd

CSVS = {
    'UVES_master':                    'data_interim/qc/uves_master.csv',
    'UVES_master_terms':              'data_interim/qc/uves_master_terms.csv',
    'UVES_relative':                  'data_interim/qc/uves_relative.csv',
    'UVES_relative_decorrelated':     'data_interim/qc/uves_relative_decorrelated.csv',
    'CHIRON_v2_master':               'data_interim/qc/chiron_v2_master.csv',
    'CHIRON_v2_decorrelated':         'data_interim/qc/chiron_v2_decorrelated.csv',
    'CHIRON_master_terms_qc':         'data_interim/qc/chiron_master_terms_qc.csv',
}

for name, path in CSVS.items():
    print("="*78)
    print(f"  {name}")
    print(f"  {path}")
    print("="*78)
    if not os.path.exists(path):
        print("  [SKIP] no existe")
        continue
    df = pd.read_csv(path)
    print(f"  N filas: {len(df)}")
    print(f"  Columnas ({len(df.columns)}): {list(df.columns)}")
    
    # Buscar columnas de tiempo y RV con nombres tipicos
    time_candidates = [c for c in df.columns
                       if any(k in c.lower() for k in ['bjd','jd','mjd','time'])]
    rv_candidates = [c for c in df.columns
                     if c.lower() in ['rv','rv_ms','rv_mlc','rv_kms','rel_rv','rv_rel']]
    err_candidates = [c for c in df.columns
                      if any(k in c.lower() for k in ['rv_err','erv','rv_uncertainty'])]
    
    print(f"  Tiempo candidatos: {time_candidates}")
    print(f"  RV candidatos:     {rv_candidates}")
    print(f"  RV err candidatos: {err_candidates}")
    
    # Stats basicos
    if time_candidates:
        tcol = time_candidates[0]
        t = df[tcol].dropna()
        if len(t) > 0:
            print(f"  Rango {tcol}: {t.min():.2f} -> {t.max():.2f}")
    
    if rv_candidates:
        rvcol = rv_candidates[0]
        rv = df[rvcol].dropna()
        if len(rv) > 0:
            print(f"  {rvcol}: median={rv.median():.4f}  std={rv.std():.4f}  N_finite={len(rv)}")
    
    print()
