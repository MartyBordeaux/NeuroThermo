from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_config
from .io import make_manifest_from_qc
from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="neurothermo-phenotypes")
    sub = parser.add_subparsers(dest="command", required=True)
    manifest = sub.add_parser("make-manifest", help="Build an analysis manifest from the frozen QC workbook")
    manifest.add_argument("--qc-xlsx", required=True)
    manifest.add_argument("--raw-root", required=True)
    manifest.add_argument("--cohort-csv")
    manifest.add_argument("--expected-included", type=int)
    manifest.add_argument("--out", required=True)
    run = sub.add_parser("run", help="Run the complete model-free phenotype analysis")
    run.add_argument("--manifest", required=True)
    run.add_argument("--config", required=False)
    run.add_argument("--data-root", required=False, help="Override input.data_root from the configuration")
    run.add_argument("--curated-events", required=False)
    run.add_argument("--curated-sweeps", required=False)
    run.add_argument("--peak-overrides", required=False)
    run.add_argument("--threshold-brackets", required=False)
    run.add_argument("--out", required=True)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "make-manifest":
        frame = make_manifest_from_qc(
            args.qc_xlsx, args.raw_root,
            cohort_csv=args.cohort_csv,
            expected_included=args.expected_included,
        )
        out = Path(args.out).expanduser()
        out.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(out, index=False)
        print(f"Wrote {len(frame)} cells to {out}; included={int(frame.include.sum())}")
    elif args.command == "run":
        cfg = load_config(args.config)
        if args.data_root:
            cfg["input"]["data_root"] = args.data_root
        overrides = {
            "curated_events_csv": args.curated_events,
            "curated_sweeps_manifest_csv": args.curated_sweeps,
            "curated_peak_overrides_csv": args.peak_overrides,
            "curated_threshold_brackets_csv": args.threshold_brackets,
        }
        for key, value in overrides.items():
            if value:
                cfg["input"][key] = value
        result = run_pipeline(args.manifest, args.out, cfg)
        print(f"Analysis completed: {result}")
