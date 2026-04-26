"""
Verifica que los FITS de ESPRESSO descargados contienen los datos
que SERVAL necesita: S2D (espectro 2D por orden) o S1D.

SERVAL clásico necesita S2D. Si solo tenemos S1D, hay que decidir:
  - usar SERVAL en modo S1D (limitado, no todos los indicadores)
  - o re-descargar productos S2D del archivo ESO
"""
from pathlib import Path
from astropy.io import fits
import pandas as pd

WORKSPACE = Path("data_work/serval/espresso_s1d")
FITS_DIR = WORKSPACE / "fits"
files = sorted(FITS_DIR.glob("*.fits"))
print(f"Found {len(files)} FITS in workspace")

# Inspeccionar 5 archivos (variados): primero, último, y 3 del medio
indices_to_check = [0, len(files)//4, len(files)//2, 3*len(files)//4, len(files)-1]
indices_to_check = sorted(set(indices_to_check))

print("\n" + "="*72)
print("INSPECCION ESTRUCTURA FITS")
print("="*72)

structure_summary = {}

for idx in indices_to_check:
    f = files[idx]
    print(f"\n[{idx}] {f.name}")
    
    try:
        with fits.open(f) as hdul:
            print(f"  N_HDUs: {len(hdul)}")
            for i_hdu, hdu in enumerate(hdul):
                hdu_type = type(hdu).__name__
                hdu_name = hdu.name if hdu.name else "(no name)"
                shape = hdu.data.shape if hdu.data is not None else "no data"
                print(f"    HDU{i_hdu}: type={hdu_type}, name='{hdu_name}', shape={shape}")
            
            # Header del HDU principal
            hdr = hdul[0].header
            
            # Buscar palabras clave de instrumento y modo
            keys_of_interest = ['INSTRUME', 'OBJECT', 'DPR.TYPE', 'PRO.CATG',
                                'PRO.TECH', 'TPL.ID', 'INS.MODE',
                                'ESO INS MODE', 'ESO PRO CATG', 'ESO PRO TECH']
            for k in keys_of_interest:
                # Astropy normaliza HIERARCH ESO claves a 'ESO INS MODE' formato
                v = hdr.get(k, hdr.get(k.replace('ESO ', ''), 'N/A'))
                if v != 'N/A':
                    print(f"    {k}: {v}")
            
            # Determinar formato (S1D o S2D)
            # Heurística: S2D tiene ~170 órdenes (ESPRESSO HR), S1D tiene 1 espectro largo
            primary_data = hdul[1].data if len(hdul) > 1 else None
            if primary_data is not None and hasattr(primary_data, 'shape'):
                if len(primary_data.shape) == 2 and primary_data.shape[0] > 50:
                    fmt = "S2D (echelle 2D)"
                elif len(primary_data.shape) == 1:
                    fmt = "S1D (1D combined)"
                else:
                    fmt = f"unknown shape {primary_data.shape}"
                print(f"  → Formato: {fmt}")
                structure_summary[f.name] = fmt
            elif primary_data is None:
                # Buscar en otros HDUs
                for i_hdu in range(1, len(hdul)):
                    if hdul[i_hdu].data is not None:
                        d = hdul[i_hdu].data
                        if hasattr(d, 'names'):  # FITS table
                            print(f"  HDU{i_hdu} es FITS table con cols: {d.names[:8]}")
                            structure_summary[f.name] = f"TABLE_{i_hdu}"
                        break
    
    except Exception as e:
        print(f"  ERROR leyendo: {e}")

print("\n" + "="*72)
print("RESUMEN DE FORMATOS DETECTADOS")
print("="*72)
from collections import Counter
fmt_count = Counter(structure_summary.values())
for fmt, n in fmt_count.most_common():
    print(f"  {fmt}: {n} archivos (de {len(structure_summary)} inspeccionados)")
