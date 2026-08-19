import argparse
import json

from . import __version__
from .config import load_config, resolve_output, resolve_upstream
from .pipeline import run_pipeline, validate_inputs


def build_parser():
    parser = argparse.ArgumentParser(
        prog="neurothermo-dependency-aware",
        description="Dependency-aware cell-level NeuroThermo robustness audit; each recorded cell is independent.",
    )
    parser.add_argument("--config", required=True, help="YAML configuration file")
    parser.add_argument("--upstream-dir", help="Override v0.3.1 output directory")
    parser.add_argument("--output-dir", help="Override results directory")
    parser.add_argument("--validate-only", action="store_true", help="Validate inputs without analysis")
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    config = load_config(args.config)
    upstream = resolve_upstream(config, args.upstream_dir)
    output = resolve_output(config, args.output_dir)
    if args.validate_only:
        print(json.dumps(validate_inputs(config, upstream), indent=2, ensure_ascii=False))
        return 0
    run_pipeline(config, upstream, output)
    print("Results: {}".format(output))
    return 0
