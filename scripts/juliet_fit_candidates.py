#!/usr/bin/env python3
"""
Juliet candidate fits for Alpha Cen B.

Models:
  M0     : offsets + jitters only
  C53    : one Keplerian around 45-60 d
  C70    : one Keplerian around 60-80 d
  C280   : one Keplerian around 220-340 d
  C53C70 : two Keplerians, C53 + C70

Default input:
  data_interim/qc/multinstrument_rvonly.csv

Units:
  Input RVs are assumed to be km/s by default and converted to m/s after
  subtracting the per-instrument median.

Example:
  python scripts/juliet_fit_candidates.py C53 --dataset all
  python scripts/juliet_fit_candidates.py C70 --dataset espresso
  python scripts/juliet_fit_candidates.py C280 --dataset all
"""

import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
import juliet


INPUT_DEFAULT = "data_interim/qc/multinstrument_rvonly.csv"


MODEL_DEFS = {
    "M0": [],
    "C53": [
        {"name": "p1", "period": [45.0, 60.0], "label": "C53"}
    ],
    "C70": [
        {"name": "p1", "period": [60.0, 80.0], "label": "C70"}
    ],
    "C280": [
        {"name": "p1", "period": [220.0, 340.0], "label": "C280_0p8AU"}
    ],
    "C53C70": [
        {"name": "p1", "period": [45.0, 60.0], "label": "C53"},
        {"name": "p2", "period": [60.0, 80.0], "label": "C70"},
    ],
}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("model", choices=sorted(MODEL_DEFS.keys()))
    ap.add_argument("--input", default=INPUT_DEFAULT)
    ap.add_argument("--dataset", choices=["all", "espresso", "uves"], default="all")
    ap.add_argument("--outbase", default="juliet_candidate_runs")
    ap.add_argument("--rv-unit", choices=["kmps", "mps"], default="kmps")
    ap.add_argument("--nlive", type=int, default=500)
    ap.add_argument("--sampler", default="dynesty")
    ap.add_argument("--force", action="store_true", help="Allow Juliet to overwrite/reuse folder.")
    return ap.parse_args()


def load_data(path, dataset, rv_unit):
    df = pd.read_csv(path)

    required = ["instrument", "bjd_tdb", "rv", "rv_err"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {missing}. Available: {list(df.columns)}")

    if "quality_flag" in df.columns:
        df = df[df["quality_flag"].fillna(0).astype(int) == 0].copy()

    df = df.dropna(subset=required).copy()

    if dataset == "espresso":
        df = df[df["instrument"].str.upper() == "ESPRESSO"].copy()
    elif dataset == "uves":
        df = df[df["instrument"].str.upper() == "UVES"].copy()

    if len(df) < 5:
        raise ValueError(f"Too few points after filtering: N={len(df)}")

    # Center per instrument before passing to Juliet.
    df["rv_centered"] = df["rv"] - df.groupby("instrument")["rv"].transform("median")

    if rv_unit == "kmps":
        df["rv_model"] = df["rv_centered"] * 1000.0
        df["rv_err_model"] = df["rv_err"] * 1000.0
    else:
        df["rv_model"] = df["rv_centered"]
        df["rv_err_model"] = df["rv_err"]

    # Sanity: avoid zero or negative errors.
    mederr = np.nanmedian(df.loc[df["rv_err_model"] > 0, "rv_err_model"])
    df["rv_err_model"] = np.where(
        df["rv_err_model"] > 0,
        df["rv_err_model"],
        mederr
    )

    return df


def make_rv_dicts(df):
    t_rv = {}
    y_rv = {}
    yerr_rv = {}

    for inst, sub in df.groupby("instrument"):
        key = str(inst)
        t_rv[key] = sub["bjd_tdb"].to_numpy(float)
        y_rv[key] = sub["rv_model"].to_numpy(float)
        yerr_rv[key] = sub["rv_err_model"].to_numpy(float)

    return t_rv, y_rv, yerr_rv


def add_instrument_priors(priors, instruments):
    for inst in instruments:
        # RVs are already median-centered. Keep offsets flexible but not huge.
        priors[f"mu_{inst}"] = {
            "distribution": "normal",
            "hyperparameters": [0.0, 10.0],
        }
        priors[f"sigma_w_{inst}"] = {
            "distribution": "loguniform",
            "hyperparameters": [1e-3, 20.0],
        }


def add_planet_priors(priors, planet_name, period_range, tref):
    pmin, pmax = period_range

    # Period.
    priors[f"P_{planet_name}"] = {
        "distribution": "uniform",
        "hyperparameters": [pmin, pmax],
    }

    # T0 over one broad representative cycle near tref.
    # For RV-only Keplerians this is mainly phase.
    priors[f"t0_{planet_name}"] = {
        "distribution": "uniform",
        "hyperparameters": [tref, tref + pmax],
    }

    # RV semi-amplitude in m/s.
    priors[f"K_{planet_name}"] = {
        "distribution": "uniform",
        "hyperparameters": [0.0, 10.0],
    }

    # Conservative eccentricity prior.
    priors[f"ecc_{planet_name}"] = {
        "distribution": "uniform",
        "hyperparameters": [0.0, 0.6],
    }

    priors[f"omega_{planet_name}"] = {
        "distribution": "uniform",
        "hyperparameters": [0.0, 360.0],
    }


def build_priors(df, model):
    priors = {}
    instruments = sorted(df["instrument"].astype(str).unique())
    add_instrument_priors(priors, instruments)

    tref = float(np.nanmin(df["bjd_tdb"]))

    for planet in MODEL_DEFS[model]:
        add_planet_priors(
            priors,
            planet_name=planet["name"],
            period_range=planet["period"],
            tref=tref,
        )

    return priors


def summarize_results(results, model, dataset, out_folder):
    post = results.posteriors
    lnz = post.get("lnZ", np.nan)
    lnzerr = post.get("lnZerr", np.nan)

    samples = post.get("posterior_samples", {})

    rows = {
        "model": model,
        "dataset": dataset,
        "lnZ": lnz,
        "lnZerr": lnzerr,
        "out_folder": str(out_folder),
    }

    for key, values in samples.items():
        if key.startswith(("P_", "K_", "ecc_", "omega_", "t0_")):
            arr = np.asarray(values, dtype=float)
            rows[f"{key}_median"] = np.nanmedian(arr)
            rows[f"{key}_std"] = np.nanstd(arr)
            rows[f"{key}_p16"] = np.nanpercentile(arr, 16)
            rows[f"{key}_p84"] = np.nanpercentile(arr, 84)

    return rows


def main():
    args = parse_args()

    df = load_data(args.input, args.dataset, args.rv_unit)

    out_folder = Path(args.outbase) / f"{args.dataset}_{args.model}"
    out_folder.mkdir(parents=True, exist_ok=True)

    t_rv, y_rv, yerr_rv = make_rv_dicts(df)
    priors = build_priors(df, args.model)

    # Save reproducibility inputs.
    df.to_csv(out_folder / "input_rv_used.csv", index=False)
    with open(out_folder / "priors_used.json", "w") as f:
        json.dump(priors, f, indent=2)

    print("\n=== DATA ===")
    print("Model:", args.model)
    print("Dataset:", args.dataset)
    print("N total:", len(df))
    print(df["instrument"].value_counts().to_string())
    print("Output:", out_folder)

    print("\n=== PRIORS ===")
    for k, v in priors.items():
        print(k, v)

    dataset = juliet.load(
        priors=priors,
        t_rv=t_rv,
        y_rv=y_rv,
        yerr_rv=yerr_rv,
        out_folder=str(out_folder),
    )

    # Juliet accepts dynesty kwargs here.
    results = dataset.fit(
        sampler=args.sampler,
        n_live_points=args.nlive,
    )

    summary = summarize_results(results, args.model, args.dataset, out_folder)

    summary_path = out_folder / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # Append global CSV.
    global_csv = Path(args.outbase) / "candidate_evidence_summary.csv"
    row = pd.DataFrame([summary])

    if global_csv.exists():
        old = pd.read_csv(global_csv)
        new = pd.concat([old, row], ignore_index=True)
        # Drop exact duplicate model/dataset rows keeping last.
        new = new.drop_duplicates(subset=["model", "dataset", "out_folder"], keep="last")
    else:
        new = row

    new.to_csv(global_csv, index=False)

    print("\n=== RESULT ===")
    print("MODEL:", args.model)
    print("DATASET:", args.dataset)
    print("lnZ:", summary["lnZ"])
    print("lnZerr:", summary["lnZerr"])

    for k, v in summary.items():
        if k.startswith("P_") or k.startswith("K_"):
            print(k, v)

    print("\nSaved:")
    print(summary_path)
    print(global_csv)


if __name__ == "__main__":
    main()
