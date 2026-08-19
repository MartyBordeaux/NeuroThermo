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


def resolve_upstream(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    requested = override if override is not None else config.get("upstream_dir", "auto")
    package_root = Path(config["_config_path"]).parent.parent
    if requested != "auto":
        candidate = Path(str(requested)).expanduser().resolve()
        if not _has_files(candidate, REQUIRED_UPSTREAM):
            missing = [str(candidate / item) for item in REQUIRED_UPSTREAM if not (candidate / item).is_file()]
            raise FileNotFoundError("Required v0.3.1 inputs are missing:\n" + "\n".join(missing))
        return candidate
    home = Path.home()
    candidates = [
        package_root.parent / "NeuroThermo_thermodynamic_phenotypes_v0_3_1" / "results" / "thermodynamic_phenotypes",
        package_root.parent / "results_v0_3_1" / "results" / "thermodynamic_phenotypes",
        home / "neurothermo" / "THERMO" / "NeuroThermo_thermodynamic_phenotypes_v0_3_1" / "results" / "thermodynamic_phenotypes",
        home / "neurothermo" / "THERMO" / "results_v0_3_1" / "results" / "thermodynamic_phenotypes",
    ]
    for candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if _has_files(candidate, REQUIRED_UPSTREAM):
            return candidate
    raise FileNotFoundError(
        "Could not auto-discover v0.3.1 results. Checked:\n" +
        "\n".join(str(item.expanduser().resolve()) for item in candidates)
    )


def resolve_raw_root(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    requested = override if override is not None else config.get("raw_root", "~/neurothermo")
    return Path(str(requested)).expanduser().resolve()


def resolve_output(config: Dict[str, Any], override: Optional[str] = None) -> Path:
    requested = override if override is not None else config.get(
        "output_dir", "results/nonoverlap_iaaft"
    )
    path = Path(str(requested)).expanduser()
    if not path.is_absolute():
        path = Path(config["_config_path"]).parent.parent / path
    return path.resolve()
