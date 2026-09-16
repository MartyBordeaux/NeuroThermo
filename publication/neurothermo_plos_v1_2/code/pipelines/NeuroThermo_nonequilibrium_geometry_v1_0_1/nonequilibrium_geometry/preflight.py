from __future__ import annotations

import numpy as np
import pandas as pd

from .data import scenario_arrays
from .model import deterministic_trace, fixed_points


def run_preflight(scenarios, cfg):
    fixed_rows, continuation_rows, endpoint_rows, membership_rows = [], [], [], []
    p_endpoints = np.asarray([0.0, 1.0])
    continuation_count = int(cfg["preflight"]["current_points"])
    current_scale = np.linspace(float(cfg["preflight"]["current_scale_min"]), float(cfg["preflight"]["current_scale_max"]), continuation_count)
    unique_endpoints = {}
    for _, row in scenarios.iterrows():
        theta_end, current_end = scenario_arrays(row, p_endpoints)
        for endpoint, theta, current in zip(("WT", "SCA3"), theta_end, current_end):
            key = (endpoint,) + tuple(np.round(np.r_[theta, current], 14))
            if key not in unique_endpoints:
                unique_endpoints[key] = {
                    "endpoint": endpoint, "theta": theta, "current": current,
                    "scenario_ids": [int(row.scenario_id)],
                }
            else:
                unique_endpoints[key]["scenario_ids"].append(int(row.scenario_id))
    endpoint_counters = {"WT": 0, "SCA3": 0}
    for item in unique_endpoints.values():
        endpoint, theta, current = item["endpoint"], item["theta"], item["current"]
        endpoint_group_id = f"{endpoint}_E{endpoint_counters[endpoint]:03d}"
        endpoint_counters[endpoint] += 1
        scenario_id = min(item["scenario_ids"])
        shared = len(item["scenario_ids"])
        for member_scenario_id in sorted(item["scenario_ids"]):
            pair_key = str(
                scenarios.loc[scenarios["scenario_id"].astype(int).eq(member_scenario_id), "biological_pair_key"].iloc[0]
            )
            membership_rows.append({
                "endpoint_group_id": endpoint_group_id,
                "endpoint": endpoint,
                "representative_scenario_id": scenario_id,
                "scenario_id": member_scenario_id,
                "biological_pair_key": pair_key,
            })
        roots = fixed_points(theta, current, cfg)
        for root_index, (x, y, z, eig) in enumerate(roots):
            fixed_rows.append({
                    "endpoint_group_id": endpoint_group_id,
                    "representative_scenario_id": scenario_id, "n_scenarios_sharing_endpoint": shared,
                    "endpoint": endpoint, "root_index": root_index,
                    "x": x, "y": y, "z": z, "max_real_eigenvalue": float(np.max(eig.real)),
                    "stable": bool(np.max(eig.real) < 0),
                    "eigenvalues": ";".join(f"{value.real:.10g}{value.imag:+.10g}j" for value in eig),
            })
        trace, ok = deterministic_trace(theta, current, cfg)
        x_range = float(np.ptp(trace[:, 0])) if len(trace) else np.nan
        endpoint_rows.append({
                "endpoint_group_id": endpoint_group_id,
                "representative_scenario_id": scenario_id, "n_scenarios_sharing_endpoint": shared,
                "endpoint": endpoint, "J": float(current),
                "finite": bool(ok), "x_range": x_range,
                "oscillatory": bool(ok and x_range >= float(cfg["preflight"]["oscillation_x_range_min"])),
                "n_fixed_points": len(roots), "n_stable_fixed_points": int(sum(np.max(item[3].real) < 0 for item in roots)),
        })
        for scale in current_scale:
            trial_current = float(current * scale)
            trial, finite = deterministic_trace(theta, trial_current, cfg)
            x_span = float(np.ptp(trial[:, 0])) if len(trial) else np.nan
            trial_roots = fixed_points(theta, trial_current, cfg)
            continuation_rows.append({
                    "endpoint_group_id": endpoint_group_id,
                    "representative_scenario_id": scenario_id, "n_scenarios_sharing_endpoint": shared,
                    "endpoint": endpoint, "J_scale": float(scale),
                    "J": trial_current, "finite": bool(finite), "x_range": x_span,
                    "oscillatory": bool(finite and x_span >= float(cfg["preflight"]["oscillation_x_range_min"])),
                    "n_fixed_points": len(trial_roots),
                    "max_stability": float(min((np.max(item[3].real) for item in trial_roots), default=np.nan)),
            })
    return (
        pd.DataFrame(fixed_rows),
        pd.DataFrame(continuation_rows),
        pd.DataFrame(endpoint_rows),
        pd.DataFrame(membership_rows),
    )
