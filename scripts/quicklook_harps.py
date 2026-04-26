import pandas as pd
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

df = pd.read_csv("data_interim/qc/harps_full_master.csv")
df = df.dropna(subset=["bjd","rv","rv_err"]).copy()

print("N =", len(df))
print(df[["bjd","rv","rv_err","berv","fwhm","contrast","snr"]].describe())

# RV vs time
plt.figure(figsize=(8,4))
plt.errorbar(df["bjd"], df["rv"], yerr=df["rv_err"], fmt='o')
plt.xlabel("BJD")
plt.ylabel("RV")
plt.title("HARPS ALPCENB RVs")
plt.tight_layout()
plt.show()

# RV vs BERV
plt.figure(figsize=(5,4))
plt.scatter(df["berv"], df["rv"])
plt.xlabel("BERV")
plt.ylabel("RV")
plt.title("RV vs BERV")
plt.tight_layout()
plt.show()

# RV vs FWHM
if df["fwhm"].notna().any():
    plt.figure(figsize=(5,4))
    plt.scatter(df["fwhm"], df["rv"])
    plt.xlabel("FWHM")
    plt.ylabel("RV")
    plt.title("RV vs FWHM")
    plt.tight_layout()
    plt.show()

# Periodogram
if len(df) >= 10:
    t = df["bjd"].values
    y = df["rv"].values - df["rv"].mean()
    dy = df["rv_err"].values
    freq, power = LombScargle(t, y, dy).autopower(
        minimum_frequency=1/5000,
        maximum_frequency=1/0.5,
        samples_per_peak=10
    )
    period = 1/freq

    plt.figure(figsize=(8,4))
    plt.plot(period, power)
    plt.xscale("log")
    plt.xlabel("Period [days]")
    plt.ylabel("Power")
    plt.title("HARPS ALPCENB Lomb-Scargle")
    plt.tight_layout()
    plt.show()
