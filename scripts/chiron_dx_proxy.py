from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits

INP = Path("data_interim/qc/chiron_working.csv")
OUT = Path("data_interim/qc/chiron_dx_proxy.csv")

def load_image(fp):
    with fits.open(fp) as hdul:
        # imagen en HDU 1
        data = hdul[1].data.astype(np.float64)
    return data

def robust_profile(img):
    # recorta bordes negros/inestables
    ny, nx = img.shape
    y0, y1 = int(0.05 * ny), int(0.95 * ny)
    x0, x1 = int(0.05 * nx), int(0.95 * nx)
    sub = img[y0:y1, x0:x1]

    # quita fondo por fila
    row_bg = np.nanmedian(sub, axis=1, keepdims=True)
    sub = sub - row_bg

    # clipping de outliers muy altos
    p99 = np.nanpercentile(sub, 99.5)
    sub = np.clip(sub, 0, p99)

    # colapso espacial -> perfil espectral 1D
    prof = np.nansum(sub, axis=0)

    # normaliza
    prof = prof - np.nanmedian(prof)
    s = np.nanstd(prof)
    if not np.isfinite(s) or s == 0:
        return prof * np.nan
    return prof / s

def xcorr_shift(ref, cur):
    n = len(ref)
    good = np.isfinite(ref) & np.isfinite(cur)
    if good.sum() < n // 2:
        return np.nan, np.nan

    a = ref[good]
    b = cur[good]

    # mismo tamaño
    m = min(len(a), len(b))
    a = a[:m]
    b = b[:m]

    fa = np.fft.rfft(a)
    fb = np.fft.rfft(b)
    c = np.fft.irfft(fa * np.conj(fb), n=m)

    # centra lags
    c = np.concatenate([c[-(m//2):], c[:m//2+1]])
    lags = np.arange(-m//2, m//2 + 1)

    i = int(np.nanargmax(c))
    lag0 = lags[i]
    peak = c[i]

    # refinamiento parabólico
    if 0 < i < len(c)-1:
        y1, y2, y3 = c[i-1], c[i], c[i+1]
        denom = (y1 - 2*y2 + y3)
        if denom != 0:
            delta = 0.5 * (y1 - y3) / denom
        else:
            delta = 0.0
    else:
        delta = 0.0

    return lag0 + delta, peak

df = pd.read_csv(INP).copy()
df = df.sort_values("date_obs").reset_index(drop=True)

profiles = []
for fp in df["file_spectrum"]:
    img = load_image(fp)
    profiles.append(robust_profile(img))

# referencia = mediana de perfiles normalizados
stack = np.vstack(profiles)
ref = np.nanmedian(stack, axis=0)
ref = ref - np.nanmedian(ref)
ref /= np.nanstd(ref)

dx = []
cc_peak = []

for prof in profiles:
    s, p = xcorr_shift(ref, prof)
    dx.append(s)
    cc_peak.append(p)

df["dx_pix_proxy"] = dx
df["ccf_peak_proxy"] = cc_peak

# centrado
df["dx_pix_proxy_demean"] = df["dx_pix_proxy"] - df["dx_pix_proxy"].mean()

df.to_csv(OUT, index=False)

print(df[[
    "stem","date_obs","bjd_tdb","airmass",
    "flux_total","x_center_mass","dx_pix_proxy","dx_pix_proxy_demean","ccf_peak_proxy"
]].to_string(index=False))
print(f"\nSaved {OUT}")
