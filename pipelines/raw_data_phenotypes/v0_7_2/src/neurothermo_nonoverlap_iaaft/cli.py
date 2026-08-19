from __future__ import annotations

import argparse

from . import __version__
from .config import load_config, resolve_output, resolve_raw_root, resolve_upstream
from .pipeline import run_pipeline, validate_inputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="neurothermo-nonoverlap-iaaft",
        description="Non-overlapping ordinal PI with validated IAAFT surrogates, NeuroThermo v0.7.2.",
    )
    parser.add_argument("--config", required=True, help="YAML configuration file")
    parser.add_argument("--upstream-dir", help="v0.3.1 thermodynamic_phenotypes directory")
    parser.add_argument("--raw-root", help="raw data root containing SCA3 and WT directories")
    parser.add_argument("--output-dir", help="output directory; defaults to results/ inside this package")
    parser.add_argument("--validate-only", action="store_true", help="validate paths/cohort and stop")
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    upstream = resolve_upstream(config, args.upstream_dir)
    raw_root = resolve_raw_root(config, args.raw_root)
    output = resolve_output(config, args.output_dir)
    if args.validate_only:
        checks = validate_inputs(config, upstream, raw_root)
        print(checks.to_string(index=False))
        return 0 if bool(checks.passed.all()) else 2
    run_pipeline(config, upstream, raw_root, output)
    print("Results written to: {}".format(output))
    return 0
