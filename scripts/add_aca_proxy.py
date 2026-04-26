from pathlib import Path
import pandas as pd
import numpy as np

def contam_proxy(sep_arcsec, seeing_arcsec=None):
    sep = np.asarray(sep_arcsec, dtype=float)
    if seeing_arcsec is None:
        # proxy mínimo: contaminación cae fuerte con separación
        return 1.0 / np.maximum(sep, 1e-3)**2
    s = np.asarray(seeing_arcsec, dtype=float)
    # proxy adimensional simple
    return np.exp(-0.5 * (sep / np.maximum(s, 1e-3))**2)

if __name__ == "__main__":
    import sys
    inp = Path(sys.argv[1])
    out = Path(sys.argv[2])

    df = pd.read_csv(inp)

    if "aca_sep_arcsec" not in df.columns:
        df["aca_sep_arcsec"] = np.nan

    seeing_col = "airmass"  # placeholder si no tienes seeing
    if seeing_col in df.columns:
        df["aca_contam_proxy"] = contam_proxy(df["aca_sep_arcsec"], None)
    else:
        df["aca_contam_proxy"] = contam_proxy(df["aca_sep_arcsec"], None)

    df.to_csv(out, index=False)
    print(f"saved {out}")
