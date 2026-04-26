"""
Fase 0 — Split HARPS pre/post fibra (mayo 2015).

Crea una nueva columna `inst_split` que distingue HARPS_pre y HARPS_post,
y guarda un nuevo master + un diagnóstico visual.

NO modifica el archivo original.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# CONFIGURACIÓN — ajusta solo si tu CSV usa otros nombres de columna
# -------------------------------------------------------------------
INPUT_CSV  = "data_interim/qc/harps_master_terms.csv"
OUTPUT_CSV = "data_interim/qc/harps_master_terms_split.csv"
DIAG_PNG   = "data_interim/qc/harps_split_diagnostic.png"

TIME_COL = "bjd"
RV_COL   = "rv"
ERV_COL  = "erv"

# Fechas del cambio de fibra HARPS (BJD)
BJD_FIBER_START = 2457161.5   # 16 mayo 2015
BJD_FIBER_END   = 2457173.5   # 28 mayo 2015

# -------------------------------------------------------------------

def main():
    if not os.path.exists(INPUT_CSV):
        sys.exit(f"[ERROR] No encuentro {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)
    print(f"[INFO] Leidas {len(df)} epocas de {INPUT_CSV}")
    print(f"[INFO] Columnas: {list(df.columns)}")

    for col in (TIME_COL, RV_COL):
        if col not in df.columns:
            sys.exit(f"[ERROR] Columna '{col}' no esta. "
                     f"Edita las variables de cabecera del script.")

    t = df[TIME_COL].values

    if t.max() < 2400000:
        print("[WARN] Los tiempos parecen MJD o JD-2400000. "
              "Aplicando correccion +2400000.5 (MJD->BJD)")
        t = t + 2400000.5
        df["_bjd_full"] = t
    else:
        df["_bjd_full"] = t

    pre_mask  = df["_bjd_full"] < BJD_FIBER_START
    post_mask = df["_bjd_full"] > BJD_FIBER_END
    gap_mask  = ~(pre_mask | post_mask)

    df["inst_split"] = np.where(pre_mask, "HARPS_pre",
                       np.where(post_mask, "HARPS_post", "HARPS_GAP"))

    n_pre  = int(pre_mask.sum())
    n_post = int(post_mask.sum())
    n_gap  = int(gap_mask.sum())

    print(f"[INFO] HARPS_pre  : {n_pre} epocas")
    print(f"[INFO] HARPS_post : {n_post} epocas")
    print(f"[INFO] HARPS_GAP  : {n_gap} epocas (deberian ser 0)")

    if n_gap > 0:
        print("[WARN] Hay epocas dentro del gap de fibra.")
        print(df.loc[gap_mask, [TIME_COL, RV_COL]])

    for sub in ("HARPS_pre", "HARPS_post"):
        sub_df = df[df["inst_split"] == sub]
        if len(sub_df) == 0:
            continue
        rv = sub_df[RV_COL].values
        baseline = sub_df['_bjd_full'].max() - sub_df['_bjd_full'].min()
        print(f"[STAT] {sub}: N={len(sub_df)}, "
              f"<RV>={np.median(rv):.3f}, "
              f"std(RV)={np.std(rv):.3f}, "
              f"baseline={baseline:.1f} d")

    df_out = df.drop(columns=["_bjd_full"])
    df_out.to_csv(OUTPUT_CSV, index=False)
    print(f"[OK]  Escrito {OUTPUT_CSV}")

    fig, ax = plt.subplots(figsize=(11, 4.5))
    for sub, color in [("HARPS_pre", "tab:blue"),
                       ("HARPS_post", "tab:red"),
                       ("HARPS_GAP", "tab:gray")]:
        m = df["inst_split"] == sub
        if m.sum() == 0:
            continue
        yerr = df[ERV_COL][m] if ERV_COL in df.columns else None
        ax.errorbar(df.loc[m, "_bjd_full"], df.loc[m, RV_COL],
                    yerr=yerr, fmt=".", ms=4, alpha=0.7,
                    color=color, label=f"{sub} (N={m.sum()})")
    ax.axvspan(BJD_FIBER_START, BJD_FIBER_END, color="orange", alpha=0.2,
               label="Fiber upgrade gap")
    ax.set_xlabel("BJD")
    ax.set_ylabel(RV_COL)
    ax.set_title("HARPS split: pre vs post fiber upgrade (May 2015)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(DIAG_PNG, dpi=120)
    print(f"[OK]  Figura {DIAG_PNG}")

if __name__ == "__main__":
    main()
