from __future__ import annotations


def _pass_fraction(frame, column, threshold):
    if frame is None or len(frame) == 0 or column not in frame:
        return None
    values = frame[column].dropna()
    return float((values <= threshold).mean()) if len(values) else None


def formalism_verdict(geometry, cfg, animal_geometry=None, oscillatory_geometry=None):
    gates = cfg["gates"]
    finite = geometry.dropna(subset=["circulation_fraction", "markov_db_violation"])
    required = float(gates["required_db_pass_fraction"])
    if finite.empty:
        adequate_fraction = continuous_pass_fraction = markov_pass_fraction = combined_pass_fraction = 0.0
    else:
        adequate = finite["stationary_current_divergence_relative"] <= float(gates["max_current_divergence_relative"])
        continuous_pass = finite["circulation_fraction"] <= float(gates["max_circulation_fraction"])
        markov_pass = finite["markov_db_violation"] <= float(gates["max_markov_db_violation"])
        adequate_fraction = float(adequate.mean())
        continuous_pass_fraction = float(continuous_pass[adequate].mean()) if adequate.any() else 0.0
        markov_pass_fraction = float(markov_pass.mean())
        combined_pass_fraction = float((adequate & continuous_pass & markov_pass).mean())
    # The continuous-current estimator is retained only as a numerical QC
    # branch. A failed divergence residual cannot be converted into evidence
    # for either equilibrium or NESS. The formalism decision therefore rests
    # on the independently estimated coarse Markov time-reversal test.
    if markov_pass_fraction < required:
        stationary_formalism = "NESS"
    else:
        stationary_formalism = "MARKOV_EQUILIBRIUM_CANDIDATE"
    equilibrium_candidate = stationary_formalism == "MARKOV_EQUILIBRIUM_CANDIDATE"
    mapping = cfg["physical_mapping"]
    mapping_complete = bool(mapping.get("enabled")) and mapping.get("beta") is not None and bool(mapping.get("energy_definition")) and bool(mapping.get("work_controls"))
    if equilibrium_candidate and mapping_complete:
        jc_status = "ELIGIBLE_FOR_SEPARATE_PHYSICAL_JARZYNSKI_CROOKS_IMPLEMENTATION"
    else:
        jc_status = "BLOCKED"
    reasons = []
    if stationary_formalism == "NESS":
        reasons.append("coarse Markov detailed-balance gate failed")
    if adequate_fraction < required:
        reasons.append("continuous-current estimator excluded by stationarity-divergence QC")
    if not mapping_complete:
        reasons.append("physical beta, energy, and conjugate work mapping not supplied")
    return {
        "stationary_formalism": stationary_formalism,
        "formalism_decision_basis": "coarse_markov_time_reversal",
        "continuous_current_used_for_formalism": False,
        "continuous_current_status": "DIAGNOSTIC_VALID" if adequate_fraction >= required else "DIAGNOSTIC_INVALID",
        "continuous_current_adequate_fraction": adequate_fraction,
        "continuous_detailed_balance_pass_fraction_when_adequate": continuous_pass_fraction,
        "markov_detailed_balance_pass_fraction": markov_pass_fraction,
        "combined_pass_fraction": combined_pass_fraction,
        "required_pass_fraction": required,
        "animal_balanced_markov_pass_fraction": _pass_fraction(
            animal_geometry, "markov_db_violation", float(gates["max_markov_db_violation"])
        ),
        "oscillatory_endpoint_markov_pass_fraction": _pass_fraction(
            oscillatory_geometry, "markov_db_violation", float(gates["max_markov_db_violation"])
        ),
        "classical_jarzynski_crooks": jc_status,
        "classical_gate_reasons": reasons,
        "hatano_sasa": "APPLICABLE_AS_NESS_IDENTITY",
        "interpretation": (
            "phi=-log(rho_ss) is an equilibrium-potential candidate only after all gates pass."
            if equilibrium_candidate
            else "phi=-log(rho_ss) is not licensed as a physical Hamiltonian; retain the NESS formulation."
        ),
    }


def validate_physical_mapping(cfg):
    mapping = cfg["physical_mapping"]
    if not mapping.get("enabled"):
        return
    missing = [key for key in ("beta", "energy_definition", "work_controls") if not mapping.get(key)]
    if missing:
        raise ValueError("physical_mapping.enabled=true but missing: " + ", ".join(missing))
    controls = {str(value) for value in mapping["work_controls"]}
    interpolated = {"b", "r", "s", "kappa_I", "J"}
    if not interpolated.issubset(controls):
        absent = sorted(interpolated - controls)
        raise ValueError(
            "Classical work cannot omit morphed controls. Supply conjugate work definitions for: "
            + ", ".join(absent)
        )
