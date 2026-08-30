from __future__ import annotations

import numpy as np

from .markers import MARKER_VARIANTS


MAIN_VIEWS = ("xyz", "xy")
PRIMARY_VIEW = "xyz"
PRIMARY_VARIANT = "seed_median_curve_isotonic"


def build_verdict(cfg, scenario_qc, ensemble, pair_markers, retention):
    gates = cfg["gates"]
    dt_values = sorted(float(value) for value in cfg["convergence"]["dt_ms"])
    fine, primary, coarse = dt_values[0], float(cfg["convergence"]["primary_dt_ms"]), dt_values[-1]
    gate_rows = []

    def add(name, passed, observed, threshold, tier="supporting"):
        gate_rows.append({"gate": name, "tier": tier, "pass": bool(passed),
                          "observed": observed, "threshold": threshold})

    expected = int(cfg["validation"]["expected_scenarios"]) * len(dt_values) * len(cfg["convergence"]["seeds"])
    add("complete_scenario_dt_seed_tasks", len(retention) == expected, int(len(retention)), expected, "fatal")
    minimum_retention = float(retention.retained_min.min())
    add("minimum_grid_retention", minimum_retention >= float(gates["minimum_grid_retention"]),
        minimum_retention, gates["minimum_grid_retention"], "fatal")

    main_qc = scenario_qc[scenario_qc.view.isin(MAIN_VIEWS)]
    for view in MAIN_VIEWS:
        fraction = float(main_qc.loc[main_qc.view == view, "pass"].mean())
        add("scenario_convergence_fraction_" + view,
            fraction >= float(gates["minimum_scenario_convergence_fraction"]),
            fraction, gates["minimum_scenario_convergence_fraction"],
            "fatal" if view == PRIMARY_VIEW else "supporting")

    lookup = ensemble.set_index(["dt_ms", "view", "marker_variant", "aggregation_order"])
    for order in ("marker_first", "curve_first"):
        for view in MAIN_VIEWS:
            for variant in MARKER_VARIANTS:
                fine_row = lookup.loc[(fine, view, variant, order)]
                core = view == PRIMARY_VIEW and variant == PRIMARY_VARIANT
                tier = "fatal" if core else "supporting"
                add("fine_pair_delta_%s_%s_%s" % (order, view, variant),
                    fine_row.median_pair_delta <= float(gates["maximum_fine_median_pair_delta"]),
                    float(fine_row.median_pair_delta), gates["maximum_fine_median_pair_delta"], tier)
                add("fine_pair_fraction_%s_%s_%s" % (order, view, variant),
                    fine_row.fraction_pairs_negative >= float(gates["minimum_fine_pair_fraction_negative"]),
                    float(fine_row.fraction_pairs_negative), gates["minimum_fine_pair_fraction_negative"], tier)
                add("fine_cell_delta_%s_%s_%s" % (order, view, variant),
                    fine_row.median_cell_delta <= float(gates["maximum_fine_median_cell_delta"]),
                    float(fine_row.median_cell_delta), gates["maximum_fine_median_cell_delta"], tier)
                add("fine_cell_fraction_%s_%s_%s" % (order, view, variant),
                    fine_row.fraction_cells_negative >= float(gates["minimum_fine_cell_fraction_negative"]),
                    float(fine_row.fraction_cells_negative), gates["minimum_fine_cell_fraction_negative"], tier)
                add("fine_cell_q75_%s_%s_%s" % (order, view, variant),
                    fine_row.q75_cell_delta < 0.0, float(fine_row.q75_cell_delta), "< 0", tier)
                add("fine_animal_fraction_%s_%s_%s" % (order, view, variant),
                    fine_row.fraction_animal_pairs_negative >= float(gates["minimum_fine_animal_pair_fraction_negative"]),
                    float(fine_row.fraction_animal_pairs_negative), gates["minimum_fine_animal_pair_fraction_negative"],
                    "supporting")
                add("fine_loo_%s_%s_%s" % (order, view, variant),
                    bool(fine_row.leave_one_animal_out_all_negative),
                    bool(fine_row.leave_one_animal_out_all_negative), True, "supporting")
                primary_row = lookup.loc[(primary, view, variant, order)]
                coarse_row = lookup.loc[(coarse, view, variant, order)]
                shift_cp = abs(float(coarse_row.median_pair_delta - primary_row.median_pair_delta))
                shift_pf = abs(float(primary_row.median_pair_delta - fine_row.median_pair_delta))
                add("ensemble_shift_coarse_primary_%s_%s_%s" % (order, view, variant),
                    shift_cp <= float(gates["ensemble_max_shift_coarse_primary"]), shift_cp,
                    gates["ensemble_max_shift_coarse_primary"], tier)
                add("ensemble_shift_primary_fine_%s_%s_%s" % (order, view, variant),
                    shift_pf <= float(gates["ensemble_max_shift_primary_fine"]), shift_pf,
                    gates["ensemble_max_shift_primary_fine"], tier)

    for dt_ms in dt_values:
        for view in MAIN_VIEWS:
            for variant in MARKER_VARIANTS:
                left = lookup.loc[(dt_ms, view, variant, "marker_first")]
                right = lookup.loc[(dt_ms, view, variant, "curve_first")]
                disagreement = abs(float(left.median_pair_delta - right.median_pair_delta))
                tier = "fatal" if view == PRIMARY_VIEW and variant == PRIMARY_VARIANT else "supporting"
                add("aggregation_order_%s_%s_%s" % (_label(dt_ms), view, variant),
                    disagreement <= float(gates["maximum_aggregation_order_disagreement"]), disagreement,
                    gates["maximum_aggregation_order_disagreement"], tier)

    enforce = bool(cfg.get("verdict", {}).get("enforce", True))
    strict_pass = all(row["pass"] for row in gate_rows)
    fatal_pass = all(row["pass"] for row in gate_rows if row["tier"] == "fatal")
    failed = [row["gate"] for row in gate_rows if not row["pass"]]
    failed_fatal = [row["gate"] for row in gate_rows if row["tier"] == "fatal" and not row["pass"]]
    failed_supporting = [row["gate"] for row in gate_rows if row["tier"] == "supporting" and not row["pass"]]
    if not enforce:
        decision = "SMOKE_ONLY"
    elif strict_pass:
        decision = "KEEP_AS_MAIN_RESULT"
    elif fatal_pass:
        decision = "KEEP_AS_ENSEMBLE_RESULT_WITH_LIMITATIONS"
    else:
        decision = "REMOVE_KL_RESULT"
    return {
        "version": "1.0.1",
        "keep_kl_as_main_result": bool(enforce and fatal_pass),
        "strict_all_gates_pass": bool(enforce and strict_pass),
        "fatal_gates_pass": bool(enforce and fatal_pass),
        "decision": decision,
        "reason": "Fatal gates protect the primary full-state ensemble result; other views, marker variants, and animal-day checks are supporting sensitivity analyses.",
        "failed_gates": failed,
        "failed_fatal_gates": failed_fatal,
        "failed_supporting_gates": failed_supporting,
        "n_gates": len(gate_rows), "n_failed_gates": len(failed),
        "n_fatal_gates": sum(row["tier"] == "fatal" for row in gate_rows),
        "n_failed_fatal_gates": len(failed_fatal),
        "n_failed_supporting_gates": len(failed_supporting),
        "gates": gate_rows,
        "interpretation": "A KEEP decision supports only an ensemble ordering along the constructed path; it does not establish a universal scenario crossing, causal precedence, irreversibility, or disease-time ordering.",
    }


def _label(value):
    return ("%.8g" % float(value)).replace(".", "p")
