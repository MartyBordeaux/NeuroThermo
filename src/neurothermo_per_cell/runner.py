from __future__ import annotations

import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize

from .data import load_observations, validate_grid
from .models import SPECS, as_dict
from .objective import loss, metrics, predict
from .splits import current_level_splits


def run(cfg: dict, resume: bool = False) -> Path:
    out = Path(cfg["output_dir"]).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    for name in ("fits", "predictions", "summary"):
        (out / name).mkdir(exist_ok=True)
    df = load_observations(cfg.get("data_path"), int(cfg["capacitance_window_ms"]))
    audit = validate_grid(df)
    selected = cfg.get("cells", "all")
    if selected != "all":
        df = df[df["cell_id"].isin(selected)].copy()
        absent = set(selected).difference(df["cell_id"].unique())
        if absent:
            raise ValueError("unknown selected cells: %s" % sorted(absent))
    tasks, rows, failures = [], [], []
    for (group, cell_id), cell in df.groupby(["group", "cell_id"], sort=True):
        for model in cfg["models"]:
            for seed in cfg["seeds"]:
                for split, train_idx, test_idx in current_level_splits(
                        cell, int(cfg["cv_folds"]), bool(cfg.get("include_full_fit", True))):
                    tasks.append((group, cell_id, model, int(seed), split,
                                  cell.loc[train_idx].copy(), cell.loc[test_idx].copy()))
    started = time.time()
    workers = max(1, int(cfg.get("workers", 1)))
    work = [(task, cfg, str(out), resume) for task in tasks]
    if workers == 1:
        results = map(_execute_task, work)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        results = executor.map(_execute_task, work, chunksize=1)
    try:
        for number, result in enumerate(results, 1):
            if result["error"]:
                failures.append({"key": result["key"], "error": result["error"]})
                print("[%d/%d] %s FAILED" % (number, len(tasks), result["key"]), flush=True)
                continue
            payload = result["payload"]
            rows.extend(_result_rows(payload, result["group"], result["cell_id"],
                                     result["model"], result["seed"], result["split"]))
            print("[%d/%d] %s %s" % (number, len(tasks), result["key"],
                                      payload["status"]), flush=True)
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
    pd.DataFrame(rows).to_csv(out / "result_rows.csv", index=False)
    pd.DataFrame(failures, columns=["key", "error"]).to_csv(out / "failures.csv", index=False)
    _summaries(out, pd.DataFrame(rows))
    status = "COMPLETE" if not failures and len(rows) == 2 * len(tasks) else "INCOMPLETE"
    metadata = {
        "pipeline_version": "0.2.1", "status": status, "audit": audit,
        "tasks_expected": len(tasks), "tasks_with_rows": len(rows) // 2,
        "failures": len(failures), "runtime_seconds": time.time() - started,
        "config": cfg,
    }
    (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return out


def _execute_task(args):
    task, cfg, out_text, resume = args
    group, cell_id, model, seed, split, train, test = task
    key = "%s__%s__%s__%s__seed%d" % (model, group, cell_id, split, seed)
    out = Path(out_text)
    fit_path = out / "fits" / (key + ".json")
    try:
        if resume and fit_path.exists():
            payload = json.loads(fit_path.read_text())
        else:
            payload = _fit_one(key, model, seed, train, test, cfg, out)
            fit_path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n")
        return {"key": key, "group": group, "cell_id": cell_id, "model": model,
                "seed": seed, "split": split, "payload": payload, "error": None}
    except Exception as exc:
        return {"key": key, "group": group, "cell_id": cell_id, "model": model,
                "seed": seed, "split": split, "payload": None, "error": repr(exc)}


def _fit_one(key, model, seed, train, test, cfg, out):
    spec, history = SPECS[model], []
    t0 = time.time()

    def callback(xk, convergence=0.0):
        history.append({"generation": len(history) + 1,
                        "best_loss": float(loss(model, xk, train, cfg)),
                        "de_convergence": float(convergence)})
        return False

    result = differential_evolution(
        lambda x: loss(model, x, train, cfg), spec.bounds, seed=seed,
        maxiter=int(cfg["optimizer"]["maxiter"]),
        popsize=int(cfg["optimizer"]["popsize"]),
        tol=float(cfg["optimizer"].get("tol", 0.005)),
        atol=float(cfg["optimizer"].get("atol", 1e-6)),
        workers=1, updating="immediate", polish=False, callback=callback,
    )
    polished = minimize(lambda x: loss(model, x, train, cfg), result.x,
                        method="Nelder-Mead", bounds=spec.bounds,
                        options={"maxiter": int(cfg["optimizer"].get("polish_maxiter", 200)),
                                 "xatol": 1e-4, "fatol": 1e-6})
    candidates = [(np.asarray(result.x, dtype=float), float(result.fun))]
    if np.isfinite(polished.fun) and _within_bounds(polished.x, spec.bounds):
        candidates.append((np.asarray(polished.x, dtype=float), float(polished.fun)))
    theta, best = min(candidates, key=lambda item: item[1])
    if not _within_bounds(theta, spec.bounds):
        raise RuntimeError("optimizer returned parameters outside declared bounds")
    plateau = _plateau(history, float(cfg["optimizer"].get("plateau_fraction", 0.005)))
    converged = bool(result.success or polished.success)
    status = "CONVERGED" if converged else "NONCONVERGED"
    train_pred = predict(model, theta, train, cfg)
    test_pred = predict(model, theta, test, cfg) if len(test) else test.copy()
    train_file = "predictions/%s__train.csv" % key
    test_file = "predictions/%s__test.csv" % key
    train_pred.to_csv(out / train_file, index=False)
    test_pred.to_csv(out / test_file, index=False)
    pop = np.asarray(result.population)
    pop_losses = np.asarray(result.population_energies)
    cutoff = best * 1.02 + 1e-12
    near = pop[pop_losses <= cutoff]
    near_cv = {}
    if len(near):
        for i, name in enumerate(spec.names):
            denom = max(abs(float(np.mean(near[:, i]))), 1e-12)
            near_cv[name] = float(np.std(near[:, i]) / denom)
    return {
        "key": key, "status": status, "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message), "polish_success": bool(polished.success),
        "polish_message": str(polished.message), "plateau_detected": plateau,
        "best_loss": best, "parameters": as_dict(model, theta),
        "parameter_order": list(spec.names), "history": history,
        "population": pop.tolist(), "population_losses": pop_losses.tolist(),
        "near_optimal_population_n": int(len(near)), "near_optimal_parameter_cv": near_cv,
        "runtime_seconds": time.time() - t0,
        "train_metrics": metrics(train_pred, cfg),
        "test_metrics": metrics(test_pred, cfg) if len(test_pred) else {},
        "train_prediction_file": train_file, "test_prediction_file": test_file,
    }


def _plateau(history, fraction):
    if len(history) < 12:
        return False
    old, new = history[-11]["best_loss"], history[-1]["best_loss"]
    return (old - new) / max(abs(old), 1e-12) < fraction


def _within_bounds(theta, bounds, atol=1e-10):
    values = np.asarray(theta, dtype=float)
    limits = np.asarray(bounds, dtype=float)
    return bool(
        np.all(np.isfinite(values))
        and np.all(values >= limits[:, 0] - atol)
        and np.all(values <= limits[:, 1] + atol)
    )


def _result_rows(payload, group, cell, model, seed, split):
    base = {"model": model, "group": group, "cell_id": cell, "seed": seed,
            "split": split, "status": payload["status"],
            "runtime_seconds": payload["runtime_seconds"],
            "parameters_json": json.dumps(payload["parameters"], sort_keys=True),
            "optimizer_success": payload["optimizer_success"],
            "plateau_detected": payload["plateau_detected"]}
    rows = []
    for partition in ("train", "test"):
        metrics_ = payload.get(partition + "_metrics", {})
        row = dict(base, partition=partition)
        row.update(metrics_)
        rows.append(row)
    return rows


def _summaries(out, rows):
    if rows.empty:
        return
    test = rows[(rows.partition == "test") & rows["fi_rmse_hz"].notna()].copy()
    group_cols = ["model", "group", "split"]
    metrics_cols = ["fi_rmse_hz", "mean_rate_error_hz", "first_spike_latency_rmse_ms",
                    "balanced_recruitment_accuracy"]
    test.groupby(group_cols)[metrics_cols].agg(["mean", "median", "std", "count"]).to_csv(
        out / "summary" / "within_cell_cv_metrics.csv")
    full = rows[(rows.partition == "train") & (rows.split == "full")].copy()
    full.to_csv(out / "summary" / "full_cell_fits.csv", index=False)
    if len(full):
        records = []
        for _, row in full.iterrows():
            for name, value in json.loads(row.parameters_json).items():
                records.append({"model": row.model, "group": row.group, "cell_id": row.cell_id,
                                "seed": row.seed, "parameter": name, "value": value,
                                "status": row.status})
        parameters = pd.DataFrame(records)
        parameters.to_csv(out / "summary" / "cell_parameters_long.csv", index=False)
        stability = (parameters.groupby(["model", "group", "cell_id", "parameter"])
                     .agg(mean=("value", "mean"), std=("value", "std"),
                          minimum=("value", "min"), maximum=("value", "max"),
                          n_seeds=("seed", "nunique"),
                          converged_fraction=("status", lambda x: float(np.mean(x == "CONVERGED"))))
                     .reset_index())
        stability["cv"] = stability["std"] / stability["mean"].abs().clip(lower=1e-12)
        stability.to_csv(out / "summary" / "parameter_stability_across_seeds.csv", index=False)
    convergence = (rows.drop_duplicates(["model", "group", "cell_id", "seed", "split"])
                   .groupby(["model", "group", "split", "status"]).size()
                   .rename("n").reset_index())
    convergence.to_csv(out / "summary" / "convergence_summary.csv", index=False)
