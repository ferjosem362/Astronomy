import pandas as pd

e = pd.read_csv("data_interim/qc/espresso_relative_detrended.csv")
u = pd.read_csv("data_interim/qc/uves_relative_decorrelated.csv")
c = pd.read_csv("data_interim/qc/chiron_v2_decorrelated.csv")

de = e[["instrument","date_obs","bjd_tdb","y_rel_resid20d"]].copy()
de = de.rename(columns={"y_rel_resid20d":"y_clean"})

du = u[["instrument","date_obs","bjd_tdb","y_rel_decor"]].copy()
du = du.rename(columns={"y_rel_decor":"y_clean"})

dc = c[["instrument","date_obs","bjd_tdb","rv_resid_decor"]].copy()
dc = dc.rename(columns={"rv_resid_decor":"y_clean"})

df = pd.concat([de, du, dc], ignore_index=True)
df.to_csv("data_interim/qc/all_relative_clean.csv", index=False)

print(df.groupby("instrument").size())
print("\nSaved data_interim/qc/all_relative_clean.csv")
