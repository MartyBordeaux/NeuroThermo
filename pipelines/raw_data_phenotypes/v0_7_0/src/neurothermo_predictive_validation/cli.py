from __future__ import annotations

import argparse

from . import __version__
from .config import load_config, resolve_output, resolve_upstream, resolve_v061
from .pipeline import run_pipeline, validate_inputs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="neurothermo-predictive-validation",
        description="Cell-level mechanistic validation of predictive dynamics for NeuroThermo v0.7.0.",
    )
    parser.add_argument("--config", required=True, help="YAML configuration file")
    parser.add_argument("--upstream-dir", help="v0.3.1 thermodynamic_phenotypes directory")
    parser.add_argument("--v061-results-dir", help="v0.6.1 current_resolved_vulnerability results directory")
    parser.add_argument("--output-dir", help="output directory; defaults to results/ inside the package")
    parser.add_argument("--validate-only", action="store_true", help="validate inputs and stop")
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    upstream = resolve_upstream(config, args.upstream_dir)
    v061 = resolve_v061(config, args.v061_results_dir)
    output = resolve_output(config, args.output_dir)
    if args.validate_only:
        checks = validate_inputs(config, upstream, v061)
        print(checks.to_string(index=False))
        return 0 if bool(checks["passed"].all()) else 2
    run_pipeline(config, upstream, v061, output)
    print("Results written to: {}".format(output))
    return 0
