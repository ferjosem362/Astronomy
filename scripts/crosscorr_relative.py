import pandas as pd
import numpy as np

df = pd.read_csv("data_interim/qc/all_relative_clean.csv")
df = df[df["y_clean"].notna() & df["bjd_tdb"].notna()].copy()

pairs = [("ESPRESSO","UVES"), ("ESPRESSO","CHIRON"), ("UVES","CHIRON")]

for a, b in pairs:
    da = df[df["instrument"] == a][["bjd_tdb","y_clean"]].dropna().sort_values("bjd_tdb")
    db = df[df["instrument"] == b][["bjd_tdb","y_clean"]].dropna().sort_values("bjd_tdb")

    xs, ys = [], []
    tb = db["bjd_tdb"].values
    yb = db["y_clean"].values

    for _, r in da.iterrows():
        t = r["bjd_tdb"]
        m = np.abs(tb - t) < 1.0
        if m.sum() == 0:
            continue
        xs.append(r["y_clean"])
        ys.append(np.nanmedian(yb[m]))

    xs = np.array(xs, float)
    ys = np.array(ys, float)
    m = np.isfinite(xs) & np.isfinite(ys)

    print(f"\n=== {a} vs {b} ===")
    print("matched points:", int(m.sum()))
    if m.sum() >= 3:
        print("corr:", np.corrcoef(xs[m], ys[m])[0,1])
    else:
        print("corr: not enough overlap")
