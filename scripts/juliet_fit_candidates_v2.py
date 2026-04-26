#!/usr/bin/env python3
"""
Juliet candidate fits v2.

Usa el input limpio de juliet_runs/M0/rvs.dat:
  col0 = time
  col1 = rv clean
  col2 = rv err clean
  col3 = instrument
  col4... = linear regressors theta0..thetaN

Modelos:
  M0
  C53  : P 45-60 d
  C70  : P 65-85 d  (equivalente refinado de M74)
  C280 : P 220-340 d, ~0.8 AU
  C53C70

Ejemplos:
  python scripts/juliet_fit_candidates_v2.py M0 --dataset all
  python scripts/juliet_fit_candidates_v2.py C53 --dataset espresso
  python scripts/juliet_fit_candidates_v2.py C70 --dataset all
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import juliet


INPUT_DEFAULT = "juliet_runs/M0/rvs.dat"

MODEL_DEFS = {
    "M0": [],
    "C53": [
        ("p1", 45.0, 60.0),
    ],
    "C70": [
        ("p1", 65.0, 85.0),
    ],
    "C280": [
        ("p1", 220.0, 340.0),
    ],
    "C53C70": [
        ("p1", 45.0, 60.0),
        ("p2", 65.0, 85.0),
    ],
}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("model", choices=sorted(MODEL_DEFS))
    ap.add_argument("--input", default=INPUT_DEFAULT)
    ap.add_argument("--dataset", choices=["all", "espresso", "uves"], default="all")
    ap.add_argument("--outbase", default="juliet_candidate_runs_v2")
    ap.add_argument("--nlive", type=int, default=500)
    ap.add_argument("--sampler", default="dynesty")
    ap.add_argument("--dynamic", action="store_true", help="Use dynamic dynesty if juliet supports it")
    return ap.parse_args()


def load_clean_rvs(path, dataset):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path, sep=r"\s+", header=None, comment="#")

    if df.shape[1] < 4:
        raise ValueError(f"{path} debe tener al menos 4 columnas. Tiene {df.shape[1]}")

    base_cols = ["time", "rv", "rv_err", "instrument"]
    nreg = df.shape[1] - 4
    reg_cols = [f"theta{i}" for i in range(nreg)]
    df.columns = base_cols + reg_cols

    if dataset == "espresso":
        df = df[df["instrument"].astype(str).str.upper() == "ESPRESSO"].copy()
    elif dataset == "uves":
        df = df[df["instrument"].astype(str).str.upper() == "UVES"].copy()
    else:
        df = df[df["instrument"].astype(str).str.upper().isin(["ESPRESSO", "UVES"])].copy()

    for c in ["time", "rv", "rv_err"] + reg_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["time", "rv", "rv_err", "instrument"]).copy()
    df = df[df["rv_err"] > 0].copy()

    if len(df) < 5:
        raise ValueError(f"Demasiados pocos puntos para {dataset}: N={len(df)}")

    # Chequeo anti-basura: aquí las RV limpias deben estar en escala de pocos m/s,
    # no en miles de m/s.
    span = df["rv"].max() - df["rv"].min()
    if span > 50:
        raise RuntimeError(
            f"RV span demasiado grande ({span}). Este no parece el rvs.dat limpio."
        )

    return df, reg_cols


def make_dicts(df, reg_cols):
    t_rv = {}
    y_rv = {}
    yerr_rv = {}
    linear_regressors_rv = {}

    for inst, sub in df.groupby("instrument"):
        inst = str(inst)
        t_rv[inst] = sub["time"].to_numpy(float)
        y_rv[inst] = sub["rv"].to_numpy(float)
        yerr_rv[inst] = sub["rv_err"].to_numpy(float)

        if reg_cols:
            linear_regressors_rv[inst] = sub[reg_cols].to_numpy(float)

    return t_rv, y_rv, yerr_rv, linear_regressors_rv


def build_priors(df, reg_cols, model):
    priors = {}
    instruments = sorted(df["instrument"].astype(str).unique())

    for inst in instruments:
        priors[f"mu_{inst}"] = {
            "distribution": "uniform",
            "hyperparameters": [-5.0, 5.0],
        }
        priors[f"sigma_w_{inst}"] = {
            "distribution": "loguniform",
            "hyperparameters": [0.0001, 5.0],
        }

        for i, _ in enumerate(reg_cols):
            priors[f"theta{i}_{inst}"] = {
                "distribution": "normal",
                "hyperparameters": [0.0, 1.0],
            }

    tref = float(df["time"].min())

    for pname, pmin, pmax in MODEL_DEFS[model]:
        priors[f"P_{pname}"] = {
            "distribution": "uniform",
            "hyperparameters": [pmin, pmax],
        }
        priors[f"t0_{pname}"] = {
            "distribution": "uniform",
            "hyperparameters": [tref, tref + pmax],
        }
        priors[f"K_{pname}"] = {
            "distribution": "uniform",
            "hyperparameters": [0.0, 5.0],
        }

        # Igual que M74: circular para evitar que la excentricidad absorba ruido/window.
        priors[f"ecc_{pname}"] = {
            "distribution": "fixed",
            "hyperparameters": 0.0,
        }
        priors[f"omega_{pname}"] = {
            "distribution": "fixed",
            "hyperparameters": 90.0,
        }

    return priors


def summarize(results, model, dataset, out_folder):
    post = results.posteriors
    lnz = post.get("lnZ", np.nan)
    lnzerr = post.get("lnZerr", np.nan)
    samples = post.get("posterior_samples", {})

    row = {
        "model": model,
        "dataset": dataset,
        "lnZ": lnz,
        "lnZerr": lnzerr,
        "out_folder": str(out_folder),
    }

    for key, values in samples.items():
        if key.startswith(("P_", "K_", "t0_", "ecc_", "omega_")):
            arr = np.asarray(values, dtype=float)
            row[f"{key}_median"] = float(np.nanmedian(arr))
            row[f"{key}_std"] = float(np.nanstd(arr))
            row[f"{key}_p16"] = float(np.nanpercentile(arr, 16))
            row[f"{key}_p84"] = float(np.nanpercentile(arr, 84))

    return row


def write_prior_dat(priors, outpath):
    with open(outpath, "w") as f:
        for k, v in priors.items():
            dist = v["distribution"]
            hp = v["hyperparameters"]
            if isinstance(hp, (list, tuple)):
                hpstr = ",".join(str(x) for x in hp)
            else:
                hpstr = str(hp)
            f.write(f"{k:20s} {dist:20s} {hpstr}\n")


def main():
    args = parse_args()

    df, reg_cols = load_clean_rvs(args.input, args.dataset)

    out_folder = Path(args.outbase) / f"{args.dataset}_{args.model}"
    out_folder.mkdir(parents=True, exist_ok=True)

    t_rv, y_rv, yerr_rv, linear_regressors_rv = make_dicts(df, reg_cols)
    priors = build_priors(df, reg_cols, args.model)

    df.to_csv(out_folder / "input_clean_rvs_used.csv", index=False)
    with open(out_folder / "priors_used.json", "w") as f:
        json.dump(priors, f, indent=2)
    write_prior_dat(priors, out_folder / "priors_v2.dat")

    print("\n=== DATA ===")
    print("Model  :", args.model)
    print("Dataset:", args.dataset)
    print("Input  :", args.input)
    print("N total:", len(df))
    print("Regressors:", len(reg_cols))
    print(df["instrument"].value_counts().to_string())
    print("\nRV range:", df["rv"].min(), df["rv"].max(), "span", df["rv"].max() - df["rv"].min())
    print("ERR range:", df["rv_err"].min(), df["rv_err"].max())

    print("\n=== PRIORS ===")
    for k, v in priors.items():
        print(k, v)

    dataset = juliet.load(
        priors=priors,
        t_rv=t_rv,
        y_rv=y_rv,
        yerr_rv=yerr_rv,
        linear_regressors_rv=linear_regressors_rv,
        out_folder=str(out_folder),
    )

    fit_kwargs = {
        "sampler": args.sampler,
        "n_live_points": args.nlive,
    }

    if args.dynamic:
        fit_kwargs["dynamic"] = True

    results = dataset.fit(**fit_kwargs)

    row = summarize(results, args.model, args.dataset, out_folder)

    if not np.isfinite(row["lnZ"]):
        raise RuntimeError(f"Fit inválido: lnZ no finito: {row['lnZ']}")

    summary_path = out_folder / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(row, f, indent=2)

    global_csv = Path(args.outbase) / "candidate_evidence_summary.csv"
    newrow = pd.DataFrame([row])

    if global_csv.exists():
        old = pd.read_csv(global_csv)
        out = pd.concat([old, newrow], ignore_index=True)
        out = out.drop_duplicates(subset=["dataset", "model"], keep="last")
    else:
        out = newrow

    out.to_csv(global_csv, index=False)

    print("\n=== RESULT ===")
    print("MODEL  :", args.model)
    print("DATASET:", args.dataset)
    print("lnZ    :", row["lnZ"])
    print("lnZerr :", row["lnZerr"])

    for k, v in row.items():
        if k.startswith("P_") or k.startswith("K_"):
            print(k, v)

    print("\nSaved:")
    print(summary_path)
    print(global_csv)


if __name__ == "__main__":
    main()
