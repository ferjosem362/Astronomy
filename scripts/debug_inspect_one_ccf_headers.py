#!/usr/bin/env python3
from pathlib import Path
from astropy.io import fits
import re

files = sorted(Path("data_raw/espresso/ccf").glob("*.fits"))
if not files:
    files = sorted(Path("data_raw/espresso").glob("**/*SCCA*.fits"))
if not files:
    raise SystemExit("No encuentro CCF FITS en data_raw/espresso/ccf ni con patrón **/*SCCA*.fits")

print("N files:", len(files))

for fp in files[:3]:
    print("\n" + "=" * 120)
    print("FILE:", fp)
    with fits.open(fp, memmap=True) as hdul:
        print("N HDU:", len(hdul))
        for i, hdu in enumerate(hdul):
            print(f"\n--- HDU {i} {type(hdu).__name__} EXTNAME={hdu.header.get('EXTNAME')} ---")
            print("data shape:", getattr(hdu.data, "shape", None))
            print("columns:", getattr(getattr(hdu.data, "columns", None), "names", None))

            hits = []
            for k, v in hdu.header.items():
                ku = str(k).upper()
                if any(tok in ku for tok in ["RV", "BERV", "CCF", "BIS", "FWHM", "CONTRAST", "MASK"]):
                    hits.append((k, v))

            print("Header hits:", len(hits))
            for k, v in hits[:120]:
                print(f"{k:35s} = {v}")
