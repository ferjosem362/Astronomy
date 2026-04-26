from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits

inv = pd.read_csv("data_interim/qc/chiron_inventory.csv").copy()
files = inv["file"].astype(str).tolist()

def load_img(fp):
    with fits.open(fp) as hdul:
        return hdul[1].data.astype(float)

def robust_norm(y):
    y = y.astype(float)
    med = np.nanmedian(y)
    y = y - med
    p95 = np.nanpercentile(np.abs(y), 95)
    if np.isfinite(p95) and p95 > 0:
        y = y / p95
    return y

# imagen de referencia
imgs = [load_img(fp) for fp in files[:min(10, len(files))]]
ref = np.nanmedian(np.stack(imgs, axis=0), axis=0)

# perfil en y
prof_y = np.nanmedian(ref, axis=1)
prof_y = prof_y - np.nanmedian(prof_y)

# detector simple de máximos locales
cand = []
for i in range(2, len(prof_y)-2):
    if prof_y[i] > prof_y[i-1] and prof_y[i] > prof_y[i+1]:
        cand.append((i, prof_y[i]))

cand = sorted(cand, key=lambda x: x[1], reverse=True)

# imponer separación mínima
orders_y = []
minsep = 25
for i, v in cand:
    if all(abs(i-j) >= minsep for j in orders_y):
        orders_y.append(i)
    if len(orders_y) >= 20:
        break

orders_y = np.array(sorted(orders_y))
half_height = 6

order_specs = {}
for j, yc in enumerate(orders_y):
    specs = []
    for fp in files:
        img = load_img(fp)
        y0 = max(0, yc-half_height)
        y1 = min(img.shape[0], yc+half_height+1)
        spec = np.nansum(img[y0:y1, :], axis=0)
        specs.append(robust_norm(spec))
    order_specs[j] = np.array(specs)

templates = {j: np.nanmedian(order_specs[j], axis=0) for j in order_specs}

def shift_corr(obs, tpl, max_shift=20, oversamp=8):
    x = np.arange(len(obs), dtype=float)
    xi = np.linspace(0, len(obs)-1, len(obs)*oversamp)
    oo = np.interp(xi, x, obs)
    tt = np.interp(xi, x, tpl)

    best_s = np.nan
    best_c = -np.inf
    for s in range(-max_shift*oversamp, max_shift*oversamp+1):
        if s >= 0:
            a = oo[s:]
            b = tt[:len(tt)-s]
        else:
            a = oo[:len(oo)+s]
            b = tt[-s:]
        if len(a) < 100:
            continue
        aa = a - np.nanmean(a)
        bb = b - np.nanmean(b)
        den = np.sqrt(np.nansum(aa**2) * np.nansum(bb**2))
        if not np.isfinite(den) or den <= 0:
            continue
        c = np.nansum(aa*bb) / den
        if c > best_c:
            best_c = c
            best_s = s / oversamp
    return best_s, best_c

rows = []
for i, fp in enumerate(files):
    rec = {"file": fp, "stem": Path(fp).name.replace(".fits.fz","")}
    shifts, weights = [], []

    for j in order_specs:
        obs = order_specs[j][i]
        tpl = templates[j]
        s, c = shift_corr(obs, tpl)
        rec[f"shift_ord_{j:02d}"] = s
        rec[f"cc_ord_{j:02d}"] = c
        if np.isfinite(s) and np.isfinite(c):
            shifts.append(s)
            weights.append(max(c, 0.0)**2)

    if len(shifts) and np.sum(weights) > 0:
        shifts = np.array(shifts)
        weights = np.array(weights)
        rec["dx_combined"] = np.sum(weights * shifts) / np.sum(weights)
        rec["dx_scatter"] = np.nanstd(shifts)
        rec["n_orders_used"] = len(shifts)
    else:
        rec["dx_combined"] = np.nan
        rec["dx_scatter"] = np.nan
        rec["n_orders_used"] = 0

    rows.append(rec)
    print("done", i+1, "/", len(files), rec["stem"], rec["dx_combined"])

out = pd.DataFrame(rows)
out.to_csv("data_interim/qc/chiron_order_pseudorv_v1.csv", index=False)
print("\nSaved data_interim/qc/chiron_order_pseudorv_v1.csv")
