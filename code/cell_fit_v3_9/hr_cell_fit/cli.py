from __future__ import annotations
import argparse, json
from .config import load_config
from .pipeline import validate, run, identify_final


def main(argv=None):
    p = argparse.ArgumentParser(prog='hr-cell-fit', description='NeuroThermo threshold-constrained joint-cell fitter v3.9')
    sub = p.add_subparsers(dest='cmd', required=True)
    for name in ('validate', 'run', 'identify'):
        q = sub.add_parser(name); q.add_argument('--config', required=True)
    args = p.parse_args(argv); cfg, _ = load_config(args.config)
    if args.cmd == 'validate': out = validate(cfg)
    elif args.cmd == 'run': out = run(cfg)
    else: out = identify_final(cfg)
    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
