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
    y = np.asarray(y, dtype=float)
    med = np.nanmedian(y)
    y = y - med
    s = np.nanpercentile(np.abs(y), 95)
    if np.isfinite(s) and s > 0:
        y = y / s
    return y

def smooth1d(y, w=9):
    if w < 3:
        return y.copy()
    ker = np.ones(w, dtype=float) / w
    return np.convolve(y, ker, mode="same")

# imagen de referencia
imgs = [load_img(fp) for fp in files[:min(10, len(files))]]
ref = np.nanmedian(np.stack(imgs, axis=0), axis=0)

# perfil en y para detectar órdenes
prof_y = np.nanmedian(ref, axis=1)
prof_y = prof_y - np.nanmedian(prof_y)

cand = []
for i in range(2, len(prof_y)-2):
    if prof_y[i] > prof_y[i-1] and prof_y[i] > prof_y[i+1]:
        cand.append((i, prof_y[i]))

cand = sorted(cand, key=lambda x: x[1], reverse=True)

orders_y = []
minsep = 25
for i, v in cand:
    if all(abs(i-j) >= minsep for j in orders_y):
        orders_y.append(i)
    if len(orders_y) >= 20:
        break
orders_y = np.array(sorted(orders_y))

half_height = 4

# extraer espectros por orden
order_specs = {}
for j, yc in enumerate(orders_y):
    specs = []
    for fp in files:
        img = load_img(fp)
        y0 = max(0, yc-half_height)
        y1 = min(img.shape[0], yc+half_height+1)
        spec = np.nansum(img[y0:y1, :], axis=0)

        # quitar forma muy lenta para resaltar líneas
        low = smooth1d(spec, w=101)
        hp = spec - low
        hp = robust_norm(hp)

        specs.append(hp)
    order_specs[j] = np.array(specs)

# templates y derivadas
templates = {}
dtemplates = {}
for j in order_specs:
    tpl = np.nanmedian(order_specs[j], axis=0)
    tpl = smooth1d(tpl, w=5)
    dtpl = np.gradient(tpl)
    templates[j] = tpl
    dtemplates[j] = dtpl

rows = []
for i, fp in enumerate(files):
    rec = {"file": fp, "stem": Path(fp).name.replace(".fits.fz","")}
    shifts = []
    weights = []

    for j in order_specs:
        obs = order_specs[j][i].copy()
        tpl = templates[j]
        dtpl = dtemplates[j]

        m = np.isfinite(obs) & np.isfinite(tpl) & np.isfinite(dtpl)
        if m.sum() < 200:
            rec[f"shift_ord_{j:02d}"] = np.nan
            rec[f"w_ord_{j:02d}"] = np.nan
            continue

        r = obs[m] - tpl[m]
        g = dtpl[m]

        den = np.sum(g*g)
        if not np.isfinite(den) or den <= 0:
            dx = np.nan
            w = np.nan
        else:
            dx = np.sum(r*g) / den
            # peso simple: contenido en derivada
            w = den / len(g)

        rec[f"shift_ord_{j:02d}"] = dx
        rec[f"w_ord_{j:02d}"] = w

        if np.isfinite(dx) and np.isfinite(w) and w > 0:
            shifts.append(dx)
            weights.append(w)

    if len(shifts) > 0 and np.sum(weights) > 0:
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
out.to_csv("data_interim/qc/chiron_order_pseudorv_v2.csv", index=False)
print("\nSaved data_interim/qc/chiron_order_pseudorv_v2.csv")
