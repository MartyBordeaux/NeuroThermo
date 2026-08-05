from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

from .data import load_observations
from .models import SPECS
from .objective import loss


def run_profiles(cfg: dict) -> Path:
    source = Path(cfg["source_run"]).expanduser().resolve()
    out = Path(cfg["output_dir"]).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    df = load_observations(cfg.get("data_path"), int(cfg["capacitance_window_ms"]))
    results = []
    for fit_file in sorted((source / "fits").glob("*.json")):
        fit = json.loads(fit_file.read_text())
        parts = fit["key"].split("__")
        model, group, cell, split = parts[:4]
        if split != "full" or fit["status"] != "CONVERGED":
            continue
        if cfg.get("cells", "all") != "all" and cell not in cfg["cells"]:
            continue
        if model not in cfg["models"]:
            continue
        cell_df = df[(df.group == group) & (df.cell_id == cell)]
        spec = SPECS[model]
        optimum = np.array([fit["parameters"][n] for n in spec.names])
        for parameter_index, parameter in enumerate(spec.names):
            lo, hi = spec.bounds[parameter_index]
            grid = np.linspace(lo, hi, int(cfg["profile_points"]))
            free = [i for i in range(len(spec.names)) if i != parameter_index]
            free_bounds = [spec.bounds[i] for i in free]
            for value in grid:
                def objective(z):
                    theta = optimum.copy()
                    theta[parameter_index] = value
                    theta[free] = z
                    return loss(model, theta, cell_df, cfg)
                res = differential_evolution(objective, free_bounds,
                                             seed=int(cfg["seed"]),
                                             maxiter=int(cfg["optimizer"]["maxiter"]),
                                             popsize=int(cfg["optimizer"]["popsize"]),
                                             polish=True)
                results.append({"model": model, "group": group, "cell_id": cell,
                                "parameter": parameter, "fixed_value": float(value),
                                "profile_loss": float(res.fun), "optimizer_success": bool(res.success)})
    pd.DataFrame(results).to_csv(out / "profile_likelihood.csv", index=False)
    return out
