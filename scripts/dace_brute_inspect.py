"""
Inspector a saco bruto: imprime tipo, len y muestra de los nodos hoja.
"""
from dace_query.spectroscopy import Spectroscopy
import pandas as pd

print("[INFO] Consultando DACE...")
data = Spectroscopy.get_timeseries(
    target='HD128621',
    sorted_by_instrument=True,
    output_format='pandas'
)

print("\n" + "="*70)
print("INSPECCION BRUTA POR NIVELES")
print("="*70)

for inst, lvl1 in data.items():
    print(f"\n>>> {inst!r}  type={type(lvl1).__name__}")
    for drs, lvl2 in lvl1.items():
        print(f"    {drs!r}  type={type(lvl2).__name__}")
        for mode, lvl3 in lvl2.items():
            t = type(lvl3).__name__
            print(f"        {mode!r}  type={t}", end="")
            
            if isinstance(lvl3, pd.DataFrame):
                print(f"  N={len(lvl3)}  cols={len(lvl3.columns)}")
                if len(lvl3) > 0 and 'rjd' in lvl3.columns:
                    bjd_min = lvl3['rjd'].min() + 2400000
                    bjd_max = lvl3['rjd'].max() + 2400000
                    yr = (lvl3['rjd'].max() - lvl3['rjd'].min()) / 365.25
                    print(f"          BJD: {bjd_min:.1f} -> {bjd_max:.1f}  ({yr:.2f} anos)")
                    print(f"          RV mediana: {lvl3['rv'].median():.4f}  std: {lvl3['rv'].std():.4f}")
                    print(f"          eRV mediana: {lvl3['rv_err'].median():.6f}")
            elif isinstance(lvl3, dict) or hasattr(lvl3, 'keys'):
                # un nivel mas
                print(f"  keys={list(lvl3.keys())[:5]}")
                for k4, v4 in lvl3.items():
                    t4 = type(v4).__name__
                    print(f"            {k4!r}  type={t4}", end="")
                    if isinstance(v4, pd.DataFrame):
                        print(f"  N={len(v4)}  cols={len(v4.columns)}")
                        if len(v4) > 0 and 'rjd' in v4.columns:
                            bjd_min = v4['rjd'].min() + 2400000
                            bjd_max = v4['rjd'].max() + 2400000
                            yr = (v4['rjd'].max() - v4['rjd'].min()) / 365.25
                            print(f"              BJD: {bjd_min:.1f} -> {bjd_max:.1f}  ({yr:.2f} anos)")
                            print(f"              RV mediana: {v4['rv'].median():.4f}  std: {v4['rv'].std():.4f}")
                    elif isinstance(v4, list):
                        print(f"  list len={len(v4)}")
                        if len(v4) > 0:
                            print(f"              primer elemento type: {type(v4[0]).__name__}")
                    else:
                        print(f"  (no DF ni dict ni lista)")
            elif isinstance(lvl3, list):
                print(f"  list len={len(lvl3)}")
                if len(lvl3) > 0:
                    print(f"          primer elemento type: {type(lvl3[0]).__name__}")
                    if hasattr(lvl3[0], '__len__'):
                        print(f"          primer elemento len: {len(lvl3[0])}")
            else:
                print(f"  (nodo hoja inesperado)")
