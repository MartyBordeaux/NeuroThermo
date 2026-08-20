from __future__ import annotations
import argparse, json
from .analysis import load_config, validate, run

def main():
    ap=argparse.ArgumentParser(prog='thermo_v1_0_1')
    sub=ap.add_subparsers(dest='cmd',required=True)
    for name in ['validate','run']:
        p=sub.add_parser(name); p.add_argument('--config',required=True)
    a=ap.parse_args(); cfg=load_config(a.config)
    if a.cmd=='validate':
        info=validate(cfg); print(json.dumps(info,indent=2))
        if info['errors']: raise SystemExit(2)
    else:
        info=validate(cfg)
        if info['errors']: print(json.dumps(info,indent=2)); raise SystemExit(2)
        print(json.dumps(run(cfg),indent=2))
if __name__=='__main__': main()
