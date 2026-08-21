from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import yaml


REQUIRED_UPSTREAM = (
    "sweep_features.csv",
    "cell_scalar_phenotypes.csv",
    "cell_integrated_phenotypes.csv",
    "disease_coordinate/cell_disease_coordinate.csv",
    "analysis_manifest.json",
)

REQUIRED_V061 = (
    "cell_current_scores_target_specific.csv",
    "cell_vulnerability_summary.csv",
    "group_curve_exact_tests.csv",
    "I_exit_exact_test.csv",
    "analysis_manifest.json",
)


def load_config(path: str) -> Dict[str, Any]:
    config_path = Path(path).expanduser().resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("The YAML configuration must contain a mapping.")
    config["_config_path"] = str(config_path)
    return config


def _has_files(path: Path, required) -> bool:
    return all((path / item).is_file() for item in required)


def _resolve(config: Dict[str, Any], key: str, override: Optional[str], required, candidates):
    requested = override if override is not None else config.get(key, "auto")
    if requested != "auto":
        candidate = Path(str(requested)).expanduser().resolve()
        if not _has_files(candidate, required):
            missing = [str(candidate / item) for item in required if not (candidate / item).is_file()]
            raise FileNotFoundError("Required inputs are missing:\n" + "\n".join(missing))
        return candidate
    for candidate in candidates:
        candidate = candidate.resolve()
        if _has_files(candidate, required):
            return candidate
    raise FileNotFoundError("Could not auto-discover {}. Checked:\n{}".format(key, "\n".join(str(x.resolve()) for x in candidates)))


def resolve_upstream(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    package_root = Path(config["_config_path"]).parent.parent
    home = Path.home()
    candidates = [
        package_root.parent / "NeuroThermo_thermodynamic_phenotypes_v0_3_1" / "results" / "thermodynamic_phenotypes",
        package_root.parent / "results_v0_3_1" / "results" / "thermodynamic_phenotypes",
        home / "neurothermo" / "THERMO" / "NeuroThermo_thermodynamic_phenotypes_v0_3_1" / "results" / "thermodynamic_phenotypes",
        home / "neurothermo" / "THERMO" / "results_v0_3_1" / "results" / "thermodynamic_phenotypes",
    ]
    return _resolve(config, "upstream_dir", override, REQUIRED_UPSTREAM, candidates)


def resolve_v061(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    package_root = Path(config["_config_path"]).parent.parent
    home = Path.home()
    candidates = [
        package_root.parent / "NeuroThermo_current_resolved_vulnerability_v0_6_1" / "results" / "current_resolved_vulnerability",
        package_root.parent / "results_v0_6_1" / "results" / "current_resolved_vulnerability",
        home / "neurothermo" / "THERMO" / "NeuroThermo_current_resolved_vulnerability_v0_6_1" / "results" / "current_resolved_vulnerability",
        home / "neurothermo" / "THERMO" / "results_v0_6_1" / "results" / "current_resolved_vulnerability",
    ]
    return _resolve(config, "v061_results_dir", override, REQUIRED_V061, candidates)


def resolve_output(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    requested = override if override is not None else config.get("output_dir", "results/predictive_dynamics_validation")
    path = Path(str(requested)).expanduser()
    if not path.is_absolute():
        path = Path(config["_config_path"]).parent.parent / path
    return path.resolve()
