"""
Inspecciona en profundidad los datos DACE de Alpha Cen B.
Llega hasta los DataFrames y muestra cuentas, baseline, columnas reales.
"""
from dace_query.spectroscopy import Spectroscopy
import pandas as pd
import numpy as np

print("[INFO] Consultando DACE...")
data = Spectroscopy.get_timeseries(
    target='HD128621',
    sorted_by_instrument=True,
    output_format='pandas'
)

print("\n" + "="*70)
print("ESTRUCTURA COMPLETA")
print("="*70)

def walk(obj, depth=0, path=""):
    indent = "  " * depth
    if hasattr(obj, 'iloc') and hasattr(obj, 'columns'):
        # Es un DataFrame
        print(f"{indent}[DATAFRAME] path={path}")
        print(f"{indent}  N filas: {len(obj)}")
        print(f"{indent}  Columnas ({len(obj.columns)}):")
        for col in obj.columns:
            print(f"{indent}    {col}")
        if len(obj) > 0 and 'rjd' in obj.columns:
            bjd_min = obj['rjd'].min() + 2400000
            bjd_max = obj['rjd'].max() + 2400000
            baseline_yr = (obj['rjd'].max() - obj['rjd'].min()) / 365.25
            print(f"{indent}  Rango BJD: {bjd_min:.1f} -> {bjd_max:.1f} ({baseline_yr:.2f} anos)")
        if len(obj) > 0:
            print(f"{indent}  Primera fila:")
            for col in list(obj.columns)[:6]:
                print(f"{indent}    {col} = {obj.iloc[0][col]}")
        return
    
    if hasattr(obj, 'keys'):
        for k in obj.keys():
            new_path = f"{path}/{k}" if path else str(k)
            print(f"{indent}{k!r}")
            walk(obj[k], depth+1, new_path)

walk(data)
