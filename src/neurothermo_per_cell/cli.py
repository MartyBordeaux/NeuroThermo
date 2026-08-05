from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from .data import load_observations, validate_grid
from .profile import run_profiles
from .runner import run


def _config(path):
    with open(path, "r") as handle:
        cfg = yaml.safe_load(handle)
    base = Path(path).resolve().parent.parent
    for field in ("output_dir", "source_run"):
        if field in cfg and not Path(cfg[field]).expanduser().is_absolute():
            cfg[field] = str(base / cfg[field])
    return cfg


def main(argv=None):
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("validate", "run", "profile"):
        p = sub.add_parser(name)
        p.add_argument("--config", required=True)
        if name == "run":
            p.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    cfg = _config(args.config)
    if args.command == "validate":
        print(json.dumps(validate_grid(load_observations(cfg.get("data_path"),
                                                         cfg["capacitance_window_ms"])), indent=2))
    elif args.command == "run":
        print(run(cfg, args.resume))
    else:
        print(run_profiles(cfg))


if __name__ == "__main__":
    main()
