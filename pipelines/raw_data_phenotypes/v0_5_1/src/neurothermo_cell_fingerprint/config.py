from pathlib import Path

import yaml


REQUIRED_INPUTS = (
    "sweep_features.csv",
    "cell_scalar_phenotypes.csv",
    "cell_integrated_phenotypes.csv",
    "disease_coordinate/cell_disease_coordinate.csv",
    "analysis_manifest.json",
)


def load_config(path):
    config_path = Path(path).expanduser().resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    config["_config_path"] = str(config_path)
    config["_project_dir"] = str(config_path.parent.parent)
    return config


def _is_upstream(path):
    return path.is_dir() and all((path / item).exists() for item in REQUIRED_INPUTS)


def resolve_upstream(config, override=None):
    if override:
        path = Path(override).expanduser().resolve()
        if not _is_upstream(path):
            raise FileNotFoundError("Invalid upstream directory: {}".format(path))
        return path

    configured = config.get("upstream_dir", "auto")
    if configured and str(configured).lower() != "auto":
        path = Path(configured).expanduser()
        if not path.is_absolute():
            path = Path(config["_project_dir"]) / path
        path = path.resolve()
        if not _is_upstream(path):
            raise FileNotFoundError("Invalid upstream directory: {}".format(path))
        return path

    project = Path(config["_project_dir"])
    thermo = Path.home() / "neurothermo" / "THERMO"
    candidates = [
        thermo / "NeuroThermo_thermodynamic_phenotypes_v0_3_1" / "results" / "thermodynamic_phenotypes",
        thermo / "results_v0_3_1" / "results" / "thermodynamic_phenotypes",
        thermo / "results_v0_3_1" / "thermodynamic_phenotypes",
        project.parent / "results_v0_3_1" / "results" / "thermodynamic_phenotypes",
        project.parent / "results_v0_3_1" / "thermodynamic_phenotypes",
    ]
    for candidate in candidates:
        if _is_upstream(candidate):
            return candidate.resolve()
    raise FileNotFoundError(
        "Could not auto-discover v0.3.1 outputs. Pass --upstream-dir pointing to "
        "the directory that contains sweep_features.csv. Checked:\n" +
        "\n".join(str(x) for x in candidates)
    )


def resolve_output(config, override=None):
    raw = override if override else config.get("output_dir", "results/dependency_aware_fingerprint")
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = Path(config["_project_dir"]) / path
    return path.resolve()
