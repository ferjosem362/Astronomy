from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits

INP = Path("data_interim/qc/uves_master_raw.csv")
OUT = Path("data_interim/qc/uves_pseudorv_sobr.csv")

C_KMS = 299792.458

# región relativamente rica en líneas y lejos de bordes
WMIN = 5000.0
WMAX = 6200.0
NGRID = 12000

# búsqueda RV relativa respecto al template
RV_GRID = np.linspace(-20.0, 20.0, 801)

def read_uves_sobr(fp: Path):
    with fits.open(fp) as hdul:
        data = hdul[1].data
        wave = np.asarray(data["WAVE"][0], dtype=float)
        flux = np.asarray(data["FLUX_REDUCED"][0], dtype=float)
        err  = np.asarray(data["ERR_REDUCED"][0], dtype=float)
        qual = np.asarray(data["QUAL_REDUCED"][0], dtype=float)
        snr  = np.asarray(data["SNR_REDUCED"][0], dtype=float)

    ok = np.isfinite(wave) & np.isfinite(flux) & np.isfinite(err)
    ok &= (wave >= WMIN) & (wave <= WMAX)
    ok &= (err > 0)

    # si QUAL_REDUCED existe, preferimos píxeles no marcados
    if qual is not None:
        ok &= np.isfinite(qual)
        ok &= (qual == 0)

    wave = wave[ok]
    flux = flux[ok]
    err = err[ok]
    snr = snr[ok] if snr is not None else None

    if len(wave) < 5000:
        raise RuntimeError("demasiados pocos puntos válidos tras QC")

    return wave, flux, err, snr

def normalize_flux(w, f):
    # continuo simple con cuantíl alto + mediana
    q = np.nanpercentile(f, 95)
    if not np.isfinite(q) or q <= 0:
        q = np.nanmedian(f)
    fn = f / q

    med = np.nanmedian(fn)
    if np.isfinite(med) and med != 0:
        fn = fn / med

    # pasamos a "absorción positiva" alrededor de cero
    fn = 1.0 - fn
    return fn

def interp_to_loggrid(w, y, loglam_grid):
    logw = np.log(w)
    return np.interp(loglam_grid, logw, y, left=np.nan, right=np.nan)

def shift_template(template_flux, dpx):
    x = np.arange(len(template_flux), dtype=float)
    xp = x - dpx
    return np.interp(x, xp, template_flux, left=np.nan, right=np.nan)

def score(obs, tpl):
    ok = np.isfinite(obs) & np.isfinite(tpl)
    if ok.sum() < 3000:
        return np.nan
    o = obs[ok]
    t = tpl[ok]

    # ajuste lineal mínimo: o ~ a + b t
    A = np.vstack([np.ones_like(t), t]).T
    coef, *_ = np.linalg.lstsq(A, o, rcond=None)
    model = A @ coef
    res = o - model
    return np.nanmean(res * res)

df = pd.read_csv(INP).copy()
df = df[df["notes"].str.contains("prod_tag=SOBR", na=False)].copy().reset_index(drop=True)

loglam_grid = np.linspace(np.log(WMIN), np.log(WMAX), NGRID)
dlog = loglam_grid[1] - loglam_grid[0]

spectra = []
meta_rows = []

for _, r in df.iterrows():
    fp = Path(r["file_spectrum"])
    try:
        w, f, e, snr_arr = read_uves_sobr(fp)
        fn = normalize_flux(w, f)
        fi = interp_to_loggrid(w, fn, loglam_grid)

        if np.isfinite(fi).sum() < 5000:
            raise RuntimeError("muy pocos puntos tras interpolación")

        spectra.append(fi)
        meta_rows.append(r)
        print("OK ", fp.name)
    except Exception as exc:
        print("ERR", fp.name, exc)

if len(spectra) < 10:
    raise SystemExit("muy pocos espectros UVES válidos")

S = np.vstack(spectra)
meta = pd.DataFrame(meta_rows).reset_index(drop=True)

template = np.nanmedian(S, axis=0)

rows = []
for i in range(len(meta)):
    obs = S[i]
    vals = []

    for rv in RV_GRID:
        dloglam = np.log(1.0 + rv / C_KMS)
        dpx = dloglam / dlog
        tpl = shift_template(template, dpx)
        vals.append(score(obs, tpl))

    vals = np.asarray(vals, dtype=float)
    j = np.nanargmin(vals)
    rv_best = RV_GRID[j]

    # refinamiento parabólico local
    rv_err = np.nan
    if 1 <= j < len(RV_GRID) - 1:
        x = RV_GRID[j-1:j+2]
        y = vals[j-1:j+2]
        try:
            a, b, c = np.polyfit(x, y, 2)
            if a > 0:
                rv_best = -b / (2*a)
                rv_err = np.sqrt(1.0 / a)
        except Exception:
            pass

    rows.append({
        "instrument": "UVES",
        "stem": meta.loc[i, "stem"],
        "target_name": meta.loc[i, "target_name"],
        "date_obs": meta.loc[i, "date_obs"],
        "mjd_obs": meta.loc[i, "mjd_obs"],
        "rv_pseudo": rv_best,
        "rv_err_pseudo": rv_err,
        "snr": meta.loc[i, "snr"],
        "airmass": meta.loc[i, "airmass"],
        "file_spectrum": meta.loc[i, "file_spectrum"],
        "notes": "UVES SOBR pseudo-RV from template matching"
    })

outdf = pd.DataFrame(rows).sort_values("mjd_obs")
OUT.parent.mkdir(parents=True, exist_ok=True)
outdf.to_csv(OUT, index=False)

print(outdf.head(25).to_string(index=False))
print(f"\nSaved {OUT} rows={len(outdf)}")
