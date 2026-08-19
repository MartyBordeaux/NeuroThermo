from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import yaml

from . import __version__
from .disease_coordinate import (
    build_disease_coordinate,
    plot_disease_coordinate,
    write_disease_coordinate_report,
)
from .extract import analyse_trace
from .io import load_traces, read_manifest
from .plotting import LABELS, plot_effect_heatmap, plot_qc, plot_response_curves
from .statistics import (
    cell_scalar_phenotypes,
    compare_cell_scalars,
    compare_groups_by_current,
    compare_integrated_cells,
    compare_response_curves,
    compare_two_part_all_cells,
    current_density_support,
    integrated_cell_features,
    response_summary,
    rheobase_brackets,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _analyse(args):
    trace, cfg = args
    return analyse_trace(trace, cfg)


def _validate_frozen_event_use(features: pd.DataFrame, cfg: dict) -> None:
    input_cfg = cfg["input"]
    if not input_cfg.get("curated_events_csv"):
        return
    expected_sweeps = input_cfg.get("expected_curated_sweeps")
    if expected_sweeps is not None:
        observed = int(features["curated_frozen_sweep"].sum())
        if observed != int(expected_sweeps):
            raise ValueError(
                f"Analysed frozen curated sweep count is {observed}; expected {expected_sweeps}"
            )
    expected_events = input_cfg.get("expected_curated_spike_events")
    if expected_events is not None:
        loaded = int(features.loc[
            features["curated_frozen_sweep"], "curated_events_loaded"
        ].sum())
        if loaded != int(expected_events):
            raise ValueError(
                f"Frozen curated events attached to traces are {loaded}; expected {expected_events}"
            )
    expected_used = input_cfg.get("expected_curated_events_used")
    if expected_used is not None:
        used = int(features.loc[
            features["curated_frozen_sweep"], "curated_events_used"
        ].sum())
        if used != int(expected_used):
            raise ValueError(
                f"Frozen curated events retained inside the analysis boundary are {used}; "
                f"expected {expected_used}"
            )
    if input_cfg.get("restrict_curated_events_to_frozen_sweeps", False):
        unexpected = features[
            (~features["curated_frozen_sweep"]) & (features["n_spikes"] > 0)
        ]
        if not unexpected.empty:
            keys = unexpected[["group", "cell_id", "sweep_index"]].to_dict("records")
            raise ValueError(f"Non-frozen sweeps supplied spikes in strict production mode: {keys}")


def _finalize_event_audit(
    event_audit: pd.DataFrame, features: pd.DataFrame, cfg: dict
) -> pd.DataFrame:
    columns = [
        "event_row_id", "group", "cell_id", "sweep_index", "current_pA",
        "time_ms", "frozen_common_domain_sweep", "peak_override_action",
        "stim_start_ms", "stim_end_ms", "event_status", "used_in_analysis",
        "used_for_threshold_probe",
    ]
    if event_audit.empty:
        return pd.DataFrame(columns=columns)
    windows = features[[
        "group", "cell_id", "sweep_index", "stim_start_s", "stim_end_s",
    ]].copy()
    audit = event_audit.merge(
        windows, on=["group", "cell_id", "sweep_index"], how="left", validate="many_to_one"
    )
    audit["stim_start_ms"] = audit["stim_start_s"] * 1000.0
    audit["stim_end_ms"] = audit["stim_end_s"] * 1000.0
    tolerance_ms = float(cfg["input"].get("curated_event_boundary_tolerance_ms", 0.0))
    accepted = audit["event_status"].eq("accepted_frozen_after_override")
    inside_window = (
        audit["stim_start_ms"].notna()
        & audit["time_ms"].ge(audit["stim_start_ms"] - tolerance_ms)
        & audit["time_ms"].le(audit["stim_end_ms"] + tolerance_ms)
    )
    inside = accepted & inside_window
    audit.loc[inside, "event_status"] = "used_frozen_event"
    audit.loc[accepted & ~inside, "event_status"] = "excluded_outside_stimulus_window"
    audit["used_in_analysis"] = audit["event_status"].eq("used_frozen_event")
    audit["used_for_threshold_probe"] = inside_window & audit["event_status"].isin(
        ["used_frozen_event", "excluded_not_frozen_sweep"]
    )
    return audit[columns].sort_values(
        ["group", "cell_id", "sweep_index", "time_ms", "event_row_id"]
    ).reset_index(drop=True)


def _event_audit_counts(event_audit: pd.DataFrame) -> pd.DataFrame:
    columns = ["group", "event_status", "events"]
    if event_audit.empty:
        return pd.DataFrame(columns=columns)
    return (
        event_audit.groupby(["group", "event_status"], dropna=False)
        .size().rename("events").reset_index().sort_values(["event_status", "group"])
    )


def _strict_pass_counts(features: pd.DataFrame) -> pd.DataFrame:
    strict = features[features["qc_pass"]].copy()
    if strict.empty:
        return pd.DataFrame(columns=["group", "cells", "analysed_sweeps", "spiking", "thermo_eligible"])
    return strict.groupby("group").agg(
        cells=("cell_id", "nunique"), analysed_sweeps=("cell_id", "size"),
        spiking=("is_spiking", "sum"), thermo_eligible=("thermo_eligible", "sum"),
    ).reset_index()


def run_pipeline(manifest_path: Union[str, Path], output_dir: Union[str, Path], cfg: dict) -> Path:
    manifest_path = Path(manifest_path).expanduser()
    output_dir = Path(output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = read_manifest(manifest_path)
    _validate_cohort(manifest, cfg)
    loaded_traces, input_event_audit = load_traces(
        manifest, manifest_path, cfg, return_event_audit=True
    )
    traces, protocol_exclusions = _filter_common_protocol(loaded_traces, cfg)
    protocol_exclusions.to_csv(output_dir / "protocol_exclusions.csv", index=False)
    workers = int(cfg["runtime"].get("workers", 1))
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_analyse, [(trace, cfg) for trace in traces]))
    else:
        rows = [analyse_trace(trace, cfg) for trace in traces]
    features = pd.DataFrame(rows).sort_values(["group", "cell_id", "sweep_index"])
    _validate_frozen_event_use(features, cfg)
    features.to_csv(output_dir / "sweep_features.csv", index=False)
    event_audit = _finalize_event_audit(input_event_audit, features, cfg)
    event_audit.to_csv(output_dir / "curated_event_audit.csv", index=False)
    event_audit_counts = _event_audit_counts(event_audit)
    event_audit_counts.to_csv(output_dir / "curated_event_audit_counts.csv", index=False)

    rheobase = rheobase_brackets(features)
    rheobase.to_csv(output_dir / "rheobase_brackets.csv", index=False)
    _validate_frozen_thresholds(rheobase, cfg, manifest_path)
    scalar_cells = cell_scalar_phenotypes(manifest, rheobase)
    scalar_cells.to_csv(output_dir / "cell_scalar_phenotypes.csv", index=False)
    scalar_stats = compare_cell_scalars(scalar_cells, cfg)
    scalar_stats.to_csv(output_dir / "group_comparisons_cell_scalars.csv", index=False)
    primary_features = cfg["statistics"].get("primary_inference_features")
    diagnostic_features = cfg["statistics"].get("diagnostic_inference_features") or []
    comparisons = compare_groups_by_current(features, cfg, primary_features)
    comparisons.to_csv(output_dir / "group_comparisons_by_current_pA.csv", index=False)
    curve_stats = compare_response_curves(features, cfg, primary_features)
    curve_stats.to_csv(output_dir / "group_comparisons_response_curves.csv", index=False)
    if diagnostic_features:
        diagnostic_current = compare_groups_by_current(features, cfg, diagnostic_features)
        diagnostic_current.insert(0, "inference_role", "diagnostic")
        diagnostic_current.to_csv(output_dir / "diagnostic_group_comparisons_by_current_pA.csv", index=False)
        diagnostic_curves = compare_response_curves(features, cfg, diagnostic_features)
        diagnostic_curves.insert(0, "inference_role", "diagnostic")
        diagnostic_curves.to_csv(output_dir / "diagnostic_group_comparisons_response_curves.csv", index=False)
    integrated, integration_coverage = integrated_cell_features(features, cfg, return_coverage=True)
    integrated.to_csv(output_dir / "cell_integrated_phenotypes.csv", index=False)
    integration_coverage.to_csv(output_dir / "integration_coverage.csv", index=False)
    integrated_stats = compare_integrated_cells(integrated, cfg)
    integrated_stats.to_csv(output_dir / "group_comparisons_integrated.csv", index=False)
    two_part_stats, two_part_coverage, two_part_values = compare_two_part_all_cells(
        features, cfg
    )
    two_part_stats.to_csv(
        output_dir / "two_part_all_cell_comparisons.csv", index=False
    )
    two_part_coverage.to_csv(
        output_dir / "two_part_all_cell_coverage.csv", index=False
    )
    two_part_values.to_csv(
        output_dir / "two_part_all_cell_values.csv", index=False
    )

    coordinate_validation = pd.DataFrame()
    if cfg.get("disease_coordinate", {}).get("enabled", True):
        coordinate_dir = output_dir / "disease_coordinate"
        coordinate_dir.mkdir(parents=True, exist_ok=True)
        coordinate_scores, coordinate_reference, coordinate_validation, coordinate_features = (
            build_disease_coordinate(features, scalar_cells, cfg)
        )
        coordinate_scores.to_csv(
            coordinate_dir / "cell_disease_coordinate.csv", index=False
        )
        coordinate_reference.to_csv(
            coordinate_dir / "disease_coordinate_reference.csv", index=False
        )
        coordinate_validation.to_csv(
            coordinate_dir / "disease_coordinate_validation.csv", index=False
        )
        coordinate_features.to_csv(
            coordinate_dir / "cell_feature_scores_long.csv", index=False
        )
        coordinate_scores[[
            "group", "cell_id", "animal_id", "evidence_grade",
            "crossfit_score_mode", "crossfit_disease_burden_z",
            "crossfit_q_endpoint_unbounded",
        ]].to_csv(
            coordinate_dir / "cross_fitted_cell_scores.csv", index=False
        )
        coordinate_scores[[
            "group", "cell_id", "animal_id", "evidence_grade",
            "stability_reference_mode",
            "stability_bootstrap_iterations_requested",
            "stability_bootstrap_iterations_valid",
            "bootstrap_p_outside_WT_robust_boundary",
            "bootstrap_p_outside_observed_WT_envelope",
            "bootstrap_p_WT_exit_consensus",
        ]].to_csv(
            coordinate_dir / "bootstrap_cell_stability.csv", index=False
        )
        plot_disease_coordinate(
            coordinate_scores, coordinate_dir / "disease_coordinate.png"
        )
        write_disease_coordinate_report(
            coordinate_scores, coordinate_validation, coordinate_dir / "README.md"
        )
        sensitivity_score_rows = []
        sensitivity_validation_rows = []
        primary_capacitance = cfg["disease_coordinate"]["capacitance_feature"]
        for capacitance in cfg["disease_coordinate"].get(
            "capacitance_sensitivity_features", [primary_capacitance]
        ):
            if capacitance == primary_capacitance:
                cap_scores = coordinate_scores
                cap_validation = coordinate_validation
            else:
                sensitivity_cfg = deepcopy(cfg)
                sensitivity_cfg["disease_coordinate"]["capacitance_feature"] = capacitance
                cap_scores, _, cap_validation, _ = build_disease_coordinate(
                    features, scalar_cells, sensitivity_cfg,
                    compute_stability=False,
                )
            score_part = cap_scores[[
                "group", "cell_id", "disease_burden_z", "q_endpoint_unbounded",
                "q_endpoint_clipped_0_1", "domains_available",
                "dynamic_values_available", "dynamic_values_expected",
                "dynamic_coverage_fraction", "evidence_grade",
                "outside_observed_WT_envelope", "WT_exit_consensus_marker",
            ]].copy()
            score_part.insert(0, "capacitance_feature", capacitance)
            sensitivity_score_rows.append(score_part)
            validation_part = cap_validation.copy()
            validation_part.insert(0, "capacitance_feature", capacitance)
            sensitivity_validation_rows.append(validation_part)
        pd.concat(sensitivity_score_rows, ignore_index=True).to_csv(
            coordinate_dir / "capacitance_window_sensitivity_cell_scores.csv",
            index=False,
        )
        pd.concat(sensitivity_validation_rows, ignore_index=True).to_csv(
            coordinate_dir / "capacitance_window_sensitivity_validation.csv",
            index=False,
        )

    sensitivity_dir = output_dir / "diagnostics_strict_pass"
    sensitivity_dir.mkdir(parents=True, exist_ok=True)
    strict_features = features.copy()
    strict_features["qc_pass"] = strict_features["qc_status"].eq("PASS")
    strict_current = compare_groups_by_current(strict_features, cfg, primary_features)
    strict_current.to_csv(sensitivity_dir / "group_comparisons_by_current_pA.csv", index=False)
    strict_curves = compare_response_curves(strict_features, cfg, primary_features)
    strict_curves.to_csv(sensitivity_dir / "group_comparisons_response_curves.csv", index=False)
    strict_integrated = integrated_cell_features(strict_features, cfg)
    strict_integrated.to_csv(sensitivity_dir / "cell_integrated_phenotypes.csv", index=False)
    strict_integrated_stats = compare_integrated_cells(strict_integrated, cfg)
    strict_integrated_stats.to_csv(sensitivity_dir / "group_comparisons_integrated.csv", index=False)
    strict_two_part, strict_two_part_coverage, _ = compare_two_part_all_cells(
        strict_features, cfg
    )
    strict_two_part.to_csv(
        sensitivity_dir / "two_part_all_cell_comparisons.csv", index=False
    )
    strict_two_part_coverage.to_csv(
        sensitivity_dir / "two_part_all_cell_coverage.csv", index=False
    )
    if cfg.get("disease_coordinate", {}).get("enabled", True):
        strict_coordinate_dir = sensitivity_dir / "disease_coordinate"
        strict_coordinate_dir.mkdir(parents=True, exist_ok=True)
        strict_scores, strict_reference, strict_validation, strict_coordinate_features = (
            build_disease_coordinate(strict_features, scalar_cells, cfg)
        )
        strict_scores.to_csv(
            strict_coordinate_dir / "cell_disease_coordinate.csv", index=False
        )
        strict_reference.to_csv(
            strict_coordinate_dir / "disease_coordinate_reference.csv", index=False
        )
        strict_validation.to_csv(
            strict_coordinate_dir / "disease_coordinate_validation.csv", index=False
        )
        strict_coordinate_features.to_csv(
            strict_coordinate_dir / "cell_feature_scores_long.csv", index=False
        )
        strict_scores[[
            "group", "cell_id", "animal_id", "evidence_grade",
            "crossfit_score_mode", "crossfit_disease_burden_z",
            "crossfit_q_endpoint_unbounded",
        ]].to_csv(
            strict_coordinate_dir / "cross_fitted_cell_scores.csv", index=False
        )
        strict_scores[[
            "group", "cell_id", "animal_id", "evidence_grade",
            "stability_reference_mode",
            "stability_bootstrap_iterations_requested",
            "stability_bootstrap_iterations_valid",
            "bootstrap_p_outside_WT_robust_boundary",
            "bootstrap_p_outside_observed_WT_envelope",
            "bootstrap_p_WT_exit_consensus",
        ]].to_csv(
            strict_coordinate_dir / "bootstrap_cell_stability.csv", index=False
        )
        plot_disease_coordinate(
            strict_scores, strict_coordinate_dir / "disease_coordinate.png"
        )
        write_disease_coordinate_report(
            strict_scores, strict_validation, strict_coordinate_dir / "README.md"
        )
    _strict_pass_counts(strict_features).to_csv(sensitivity_dir / "run_counts.csv", index=False)

    curve_features = list(LABELS)
    physical_curves = response_summary(features, "current_pA", curve_features, cfg)
    physical_curves.to_csv(output_dir / "response_curves_current_pA.csv", index=False)
    j_columns = ["J_pA_per_pF", "J_10ms_pA_per_pF", "J_20ms_pA_per_pF", "J_50ms_pA_per_pF"]
    for axis in j_columns:
        summary = response_summary(features, axis, curve_features, cfg)
        summary.to_csv(output_dir / f"response_curves_{axis}.csv", index=False)
        support = current_density_support(summary)
        support.to_csv(output_dir / f"shared_support_{axis}.csv", index=False)
        if axis in {"J_pA_per_pF", "J_20ms_pA_per_pF"}:
            plot_response_curves(
                summary, output_dir / f"phenotype_curves_{axis}.png",
                "Current density, pA/pF", shared_support_only=True,
            )
    plot_response_curves(
        physical_curves, output_dir / "phenotype_curves_current_pA.png",
        "Injected current, pA", shared_support_only=True, inference=curve_stats,
    )
    plot_effect_heatmap(comparisons, output_dir / "effect_heatmap_current_pA.png")
    plot_qc(features, output_dir / "qc_summary.png")

    qc_counts = features.groupby("group").agg(
        cells=("cell_id", "nunique"), analysed_sweeps=("cell_id", "size"),
        qc_pass_clean=("qc_status", lambda x: int((x == "PASS").sum())),
        qc_warning=("qc_warning", "sum"), qc_fatal=("qc_fatal", "sum"),
        spiking=("is_spiking", "sum"),
        thermo_eligible=("thermo_eligible", "sum"),
    ).reset_index()
    excluded_counts = protocol_exclusions.groupby("group").size().rename("protocol_excluded").reset_index()
    qc_counts = qc_counts.merge(excluded_counts, on="group", how="left")
    qc_counts["protocol_excluded"] = qc_counts["protocol_excluded"].fillna(0).astype(int)
    qc_counts.to_csv(output_dir / "run_counts.csv", index=False)
    config_path = output_dir / "resolved_config.yaml"
    config_path.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
    input_paths = sorted({Path(x) for x in features.source_path.unique() if Path(x).exists()})
    auxiliary_paths = {}
    for key in [
        "curated_events_csv", "curated_sweeps_manifest_csv",
        "curated_peak_overrides_csv", "curated_threshold_brackets_csv",
        "curated_hash_manifest_json", "sweep_overrides_csv",
    ]:
        value = cfg["input"].get(key)
        if not value:
            continue
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = manifest_path.resolve().parent / path
        if path.exists():
            auxiliary_paths[key] = path
    provenance = {
        "pipeline_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": _sha256(manifest_path),
        "resolved_config_sha256": _sha256(config_path),
        "n_input_files": len(input_paths),
        "n_loaded_sweeps": len(loaded_traces),
        "n_protocol_excluded_sweeps": len(protocol_exclusions),
        "n_analysed_sweeps": len(features),
        "input_sha256": {str(p): _sha256(p) for p in input_paths},
        "auxiliary_input_sha256": {key: {"path": str(path), "sha256": _sha256(path)} for key, path in auxiliary_paths.items()},
        "analysed_spiking_sweeps": int(features["is_spiking"].sum()),
        "analysed_spike_events": int(features["n_spikes"].sum()),
        "frozen_curated_sweeps_in_common_domain": int(features["curated_frozen_sweep"].sum()),
        "frozen_curated_events_loaded": int(features.loc[features["curated_frozen_sweep"], "curated_events_loaded"].sum()),
        "frozen_curated_events_used": int(features.loc[features["curated_frozen_sweep"], "curated_events_used"].sum()),
        "frozen_curated_events_outside_stimulus_window": int(features.loc[features["curated_frozen_sweep"], "curated_events_outside_stimulus_window"].sum()),
        "curated_event_boundary_tolerance_ms": float(cfg["input"].get("curated_event_boundary_tolerance_ms", 0.0)),
        "all_cell_inference_cohort": {
            str(group): int(count)
            for group, count in manifest.groupby("group").size().items()
        },
        "disease_coordinate": {
            "enabled": bool(cfg.get("disease_coordinate", {}).get("enabled", True)),
            "domains": ["structure", "spike_timing", "predictive_dynamics"],
            "formula_frozen_from_version": "0.3.0",
            "descriptive_auc_scope": "internal_training_cohort",
            "cross_fitted_auc_scope": "leave_one_WT_out_internal_not_external",
            "stability_bootstrap_iterations": int(
                cfg.get("disease_coordinate", {}).get(
                    "stability_bootstrap_iterations",
                    cfg["statistics"]["bootstrap_iterations"],
                )
            ),
            "is_endpoint_similarity_not_time": True,
            "animal_ids_recoverable": bool(
                "animal_id" in features
                and features["animal_id"].notna().all()
                and not features["animal_id"].astype(str).eq("NA_NOT_RECOVERABLE").any()
            ),
        },
        "excluded_unfrozen_curated_events": int((event_audit["event_status"] == "excluded_not_frozen_sweep").sum()) if not event_audit.empty else 0,
        "excluded_peak_override_events": int((event_audit["event_status"] == "excluded_peak_override").sum()) if not event_audit.empty else 0,
        "interpretation_limits": [
            "External electrical work is not total metabolic energy.",
            "Voltage-path KL is an observable, coarse-grained irreversibility proxy and lower bound, not total cellular entropy production.",
            "Permutation and spectral entropies are information-theoretic descriptors, not thermodynamic entropy.",
            "Endpoint WT/SCA3 data cannot identify a biological point of no return.",
        ],
    }
    (output_dir / "analysis_manifest.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_summary(
        output_dir, qc_counts, comparisons, integrated_stats, scalar_stats,
        curve_stats, protocol_exclusions, event_audit_counts, two_part_stats,
        coordinate_validation,
    )
    return output_dir


def _validate_cohort(manifest: pd.DataFrame, cfg: dict) -> None:
    expected = cfg.get("cohort") or {}
    if not expected:
        return
    if "expected_total" in expected and len(manifest) != int(expected["expected_total"]):
        raise ValueError(f"Manifest contains {len(manifest)} included cells; expected {expected['expected_total']}")
    for group, count in (expected.get("expected_by_group") or {}).items():
        actual = int((manifest["group"] == str(group).upper()).sum())
        if actual != int(count):
            raise ValueError(f"Manifest contains {actual} {group} cells; expected {count}")


def _filter_common_protocol(traces, cfg: dict):
    input_cfg = cfg["input"]
    if not input_cfg.get("enforce_common_current_grid", False):
        return traces, pd.DataFrame(columns=[
            "group", "cell_id", "sweep_index", "current_pA", "source_path", "exclusion_reason",
        ])
    allowed = np.asarray(input_cfg.get("analysis_currents_pA") or [], float)
    if len(allowed) == 0:
        raise ValueError("enforce_common_current_grid requires analysis_currents_pA")
    tolerance = float(input_cfg.get("current_tolerance_pA", 1e-6))
    kept, excluded = [], []
    for trace in traces:
        matches = np.flatnonzero(np.abs(allowed - float(trace.current_pA)) <= tolerance)
        if len(matches) == 1:
            trace.current_pA = float(allowed[matches[0]])
            kept.append(trace)
        else:
            excluded.append({
                "group": trace.group, "cell_id": trace.cell_id,
                "sweep_index": trace.sweep_index, "current_pA": trace.current_pA,
                "source_path": trace.source_path,
                "exclusion_reason": "outside_common_current_grid" if len(matches) == 0 else "ambiguous_current_match",
            })
    audit = pd.DataFrame(excluded, columns=[
        "group", "cell_id", "sweep_index", "current_pA", "source_path", "exclusion_reason",
    ])
    if input_cfg.get("require_complete_current_grid", False):
        rows = pd.DataFrame([
            {"group": t.group, "cell_id": t.cell_id, "current_pA": t.current_pA} for t in kept
        ])
        failures = []
        for (group, cell_id), cell in rows.groupby(["group", "cell_id"]):
            observed = np.sort(cell["current_pA"].to_numpy(float))
            if len(observed) != len(allowed) or not np.allclose(observed, np.sort(allowed), atol=tolerance, rtol=0):
                failures.append((group, cell_id, observed.tolist()))
        if failures:
            raise ValueError(f"Cells lack the complete common current grid: {failures}")
    return kept, audit


def _validate_frozen_thresholds(rheobase: pd.DataFrame, cfg: dict, manifest_path: Path) -> None:
    value = cfg["input"].get("curated_threshold_brackets_csv")
    if not value:
        return
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_path.resolve().parent / path
    frozen = pd.read_csv(path)
    required = {
        "group", "cell_id", "nonspiking_current_pA", "first_spiking_current_pA",
        "threshold_bracket_width_pA",
    }
    missing = required.difference(frozen.columns)
    if missing:
        raise ValueError(f"Frozen threshold table lacks: {', '.join(sorted(missing))}")
    merged = rheobase.merge(frozen, on=["group", "cell_id"], how="outer", indicator=True)
    if not merged["_merge"].eq("both").all():
        raise ValueError("Computed and frozen threshold cohorts do not match exactly")
    checks = [
        ("rheobase_lower_nonspiking_pA", "nonspiking_current_pA"),
        ("rheobase_upper_spiking_pA", "first_spiking_current_pA"),
        ("rheobase_bracket_width_pA", "threshold_bracket_width_pA"),
    ]
    for computed, expected in checks:
        if not np.allclose(merged[computed], merged[expected], atol=1e-6, rtol=0, equal_nan=False):
            bad = merged.loc[~np.isclose(merged[computed], merged[expected], atol=1e-6, rtol=0), ["group", "cell_id", computed, expected]]
            raise ValueError(f"Computed rheobase disagrees with frozen v3.5 brackets:\n{bad.to_string(index=False)}")


def _write_summary(
    output_dir, counts, comparisons, integrated_stats, scalar_stats, curve_stats,
    protocol_exclusions, event_audit_counts, two_part_stats,
    coordinate_validation,
):
    header = "| " + " | ".join(map(str, counts.columns)) + " |"
    divider = "| " + " | ".join(["---"] * len(counts.columns)) + " |"
    body = ["| " + " | ".join(map(str, row)) + " |" for row in counts.itertuples(index=False, name=None)]
    counts_markdown = "\n".join([header, divider, *body])
    lines = [
        "# NeuroThermo model-free phenotype run",
        "",
        "## Dataset",
        "",
        counts_markdown,
        "",
        "## Primary interpretation",
        "",
        "- Current in pA is the matched experimental protocol axis.",
        "- J = I/Cm is emitted only inside shared WT/SCA3 support, including 10/20/50 ms capacitance sensitivity.",
        "- Non-spiking sweeps enter rheobase brackets; thermodynamic trace metrics require at least the configured number of spikes.",
        "- Curated Stage-1 events replace automatic peak detection in production.",
        "- Production spike events are restricted to the hash-checked frozen v3.5 sweep keys; every rejected or window-excluded event is audited.",
        "- Latency is preserved and never removed by trace alignment.",
        "- Inference is cell-level. Sweep fragments are not treated as independent biological replicates.",
        "- Integrated inference is restricted to unconditional features with complete configured common-current support.",
        "- Conditional high-current inference is two-part: all accepted cells enter binary availability, while values are compared only where physically defined.",
        "- No cell is removed globally for an incomplete conditional curve and no value is imputed on a non-spiking sweep.",
        "- The disease coordinate gives equal weight to log-capacitance, mean-ISI and predictive-information domains after robust WT-reference scaling.",
        "- q=0 and q=1 are the WT and SCA3 endpoint medians. q is not a disease probability and injected current is not disease time.",
        "- Path-KL excess is the sole primary irreversibility metric; raw and signed bias-corrected variants are diagnostics.",
        "",
        "## Scope limits",
        "",
        "External work is clamp-supplied incremental electrical work. Path KL is a bias-controlled observable irreversibility proxy. Information entropies are not heat or metabolic entropy. This endpoint analysis separates WT and SCA3 phenotypes but does not identify biological disease time or a point of no return.",
        "",
        f"Current-resolved tests: {len(comparisons)}.",
        f"Response-curve global tests: {len(curve_stats)}.",
        f"Cell-scalar tests: {len(scalar_stats)}.",
        f"Integrated cell-level tests: {len(integrated_stats)}.",
        f"All-cell two-part tests: {int(two_part_stats['permutation_p'].notna().sum()) if 'permutation_p' in two_part_stats else 0}.",
        f"Protocol-excluded sweeps: {len(protocol_exclusions)}.",
    ]
    if not coordinate_validation.empty:
        coordinate = coordinate_validation.iloc[0]
        lines.extend([
            "",
            "## Model-free disease coordinate",
            "",
            f"- Descriptive WT/SCA3 endpoint AUC: {coordinate['descriptive_auc_SCA3_vs_WT']:.4f}.",
            f"- Internal leave-one-WT-out AUC: {coordinate['internal_cross_fitted_auc_SCA3_vs_WT']:.4f}; minimum SCA3-minus-WT margin: {coordinate['crossfit_min_SCA3_minus_max_WT_margin']:.6g}.",
            f"- Exact recomputed two-sided permutation p: {coordinate['permutation_p_two_sided']:.6g}.",
            f"- SCA3 cells above the robust WT boundary: {int(coordinate['SCA3_above_WT_robust_boundary'])}/{int(coordinate['n_SCA3'])}.",
            "- Both AUC values are internal diagnostics, not external or clinical performance.",
            "- WT-exit is an operational multi-domain marker. A biological transition time requires intermediate or longitudinal samples.",
        ])
    if not event_audit_counts.empty:
        lines.extend(["", "## Curated event audit", ""])
        for row in event_audit_counts.itertuples(index=False):
            lines.append(f"- {row.group} / {row.event_status}: {int(row.events)} events.")
    (output_dir / "RUN_SUMMARY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
