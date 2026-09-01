from __future__ import annotations
import argparse, json
from pathlib import Path
import yaml
from .surface import validate, run_all


def load_cfg(path):
    p = Path(path)
    cfg = yaml.safe_load(p.read_text(encoding='utf-8'))
    base = p.parent.parent.resolve()
    data_root = Path(cfg['data']['root'])
    if not data_root.is_absolute():
        cfg['data']['root'] = str((base / data_root).resolve())
    out = Path(cfg['output']['dir'])
    if not out.is_absolute():
        cfg['output']['dir'] = str((base / out).resolve())
    return cfg


def main(argv=None):
    ap = argparse.ArgumentParser(prog='transition_v1_2')
    sp = ap.add_subparsers(dest='cmd', required=True)
    for name in ['validate', 'run']:
        p = sp.add_parser(name)
        p.add_argument('--config', required=True)
        if name == 'run':
            p.add_argument('--no-resume', action='store_true')
    ns = ap.parse_args(argv)
    cfg = load_cfg(ns.config)
    if ns.cmd == 'validate':
        print(json.dumps(validate(cfg), indent=2, sort_keys=True))
    else:
        print(json.dumps(run_all(cfg, resume=not ns.no_resume), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
