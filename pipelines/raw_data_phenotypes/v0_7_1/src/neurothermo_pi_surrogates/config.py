from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence

import yaml


REQUIRED_UPSTREAM = (
    "sweep_features.csv",
    "cell_scalar_phenotypes.csv",
    "cell_integrated_phenotypes.csv",
    "disease_coordinate/cell_disease_coordinate.csv",
    "analysis_manifest.json",
)

REQUIRED_V070 = (
    "crossfit_residualized_predictive_information.csv",
    "mode_group_curve_exact_tests.csv",
    "mode_cell_vulnerability_summary.csv",
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


def _has_files(path: Path, required: Sequence[str]) -> bool:
    return all((path / item).is_file() for item in required)


def _resolve_results(
    config: Dict[str, Any], key: str, override: Optional[str], required: Sequence[str], candidates
) -> Path:
    requested = override if override is not None else config.get(key, "auto")
    if requested != "auto":
        candidate = Path(str(requested)).expanduser().resolve()
        if not _has_files(candidate, required):
            missing = [str(candidate / item) for item in required if not (candidate / item).is_file()]
            raise FileNotFoundError("Required inputs are missing:\n" + "\n".join(missing))
        return candidate
    for candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if _has_files(candidate, required):
            return candidate
    raise FileNotFoundError(
        "Could not auto-discover {}. Checked:\n{}".format(
            key, "\n".join(str(x.expanduser().resolve()) for x in candidates)
        )
    )


def resolve_upstream(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    package_root = Path(config["_config_path"]).parent.parent
    home = Path.home()
    return _resolve_results(
        config, "upstream_dir", override, REQUIRED_UPSTREAM,
        [
            package_root.parent / "NeuroThermo_thermodynamic_phenotypes_v0_3_1" / "results" / "thermodynamic_phenotypes",
            package_root.parent / "results_v0_3_1" / "results" / "thermodynamic_phenotypes",
            home / "neurothermo" / "THERMO" / "NeuroThermo_thermodynamic_phenotypes_v0_3_1" / "results" / "thermodynamic_phenotypes",
            home / "neurothermo" / "THERMO" / "results_v0_3_1" / "results" / "thermodynamic_phenotypes",
        ],
    )


def resolve_v070(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    package_root = Path(config["_config_path"]).parent.parent
    home = Path.home()
    return _resolve_results(
        config, "v070_results_dir", override, REQUIRED_V070,
        [
            package_root.parent / "NeuroThermo_predictive_dynamics_validation_v0_7_0" / "results" / "predictive_dynamics_validation",
            package_root.parent / "results_v0_7_0" / "results" / "predictive_dynamics_validation",
            home / "neurothermo" / "THERMO" / "NeuroThermo_predictive_dynamics_validation_v0_7_0" / "results" / "predictive_dynamics_validation",
            home / "neurothermo" / "THERMO" / "results_v0_7_0" / "results" / "predictive_dynamics_validation",
        ],
    )


def resolve_raw_root(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    requested = override if override is not None else config.get("raw_root", "~/neurothermo")
    return Path(str(requested)).expanduser().resolve()


def resolve_output(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    requested = override if override is not None else config.get(
        "output_dir", "results/predictive_information_surrogates"
    )
    path = Path(str(requested)).expanduser()
    if not path.is_absolute():
        path = Path(config["_config_path"]).parent.parent / path
    return path.resolve()
