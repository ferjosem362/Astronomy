"""
Inspecciona la estructura real que devuelve DACE para Alpha Cen B.
"""
from dace_query.spectroscopy import Spectroscopy
import json

print("[INFO] Consultando DACE...")
data = Spectroscopy.get_timeseries(
    target='HD128621',
    sorted_by_instrument=True,
    output_format='pandas'
)

print(f"\n[INFO] Tipo del objeto raiz: {type(data).__name__}")
print(f"[INFO] Claves de primer nivel ({len(data)}):")
for k in data.keys():
    v = data[k]
    print(f"  {k!r} -> {type(v).__name__}")
    
    if hasattr(v, 'keys'):
        # Es otro dict/defaultdict
        for k2 in v.keys():
            v2 = v[k2]
            print(f"      {k2!r} -> {type(v2).__name__}", end="")
            if hasattr(v2, '__len__'):
                print(f" len={len(v2)}", end="")
            if hasattr(v2, 'columns'):
                print(f" cols={list(v2.columns)[:8]}...", end="")
            print()
            
            # Si es dataframe, mostramos una fila de ejemplo
            if hasattr(v2, 'iloc') and len(v2) > 0:
                row = v2.iloc[0]
                print(f"          primera fila (primeros 8 campos):")
                for col in list(v2.columns)[:8]:
                    print(f"            {col} = {row[col]}")
    elif hasattr(v, 'columns'):
        # Es directamente un DataFrame
        print(f"      cols={list(v.columns)[:8]}...")
        if len(v) > 0:
            print(f"      primera fila (primeros 8 campos):")
            for col in list(v.columns)[:8]:
                print(f"        {col} = {v.iloc[0][col]}")
