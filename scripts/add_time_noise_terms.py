from pathlib import Path
import pandas as pd
import numpy as np

YEAR = 365.25

def add_terms(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "bjd_tdb" not in out.columns:
        raise ValueError("missing bjd_tdb")
    t = out["bjd_tdb"].astype(float).values
    x = 2.0 * np.pi * t / YEAR
    out["year_sin_1"] = np.sin(x)
    out["year_cos_1"] = np.cos(x)
    out["year_sin_2"] = np.sin(2 * x)
    out["year_cos_2"] = np.cos(2 * x)
    return out

if __name__ == "__main__":
    import sys
    inp = Path(sys.argv[1])
    out = Path(sys.argv[2])
    df = pd.read_csv(inp)
    df = add_terms(df)
    df.to_csv(out, index=False)
    print(f"saved {out}")
