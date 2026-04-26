import pandas as pd
import numpy as np
from astropy.timeseries import LombScargle

df = pd.read_csv("data_interim/qc/all_relative_clean.csv")
df = df[df["y_clean"].notna() & df["bjd_tdb"].notna()].copy()

t = df["bjd_tdb"].values
y = df["y_clean"].values
e = np.ones_like(y)

ls = LombScargle(t, y, e)
freq, power = ls.autopower(
    minimum_frequency=1/3000,
    maximum_frequency=1/1.0,
    samples_per_peak=20
)

period = 1.0 / freq

for p0 in [10.9, 70.0, 74.0, 230.0, 290.0]:
    m = (period > 0.9*p0) & (period < 1.1*p0)
    d = pd.DataFrame({"period": period[m], "power": power[m]}).sort_values("power", ascending=False).head(10)
    print(f"\n=== around {p0} d ===")
    print(d.to_string(index=False))
