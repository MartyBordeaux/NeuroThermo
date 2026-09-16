from __future__ import annotations
import argparse, json
from pathlib import Path
import yaml
from .core import validate, run


def load_config(path):
    with open(path,"r",encoding="utf-8") as f:
        return yaml.safe_load(f)


def main():
    ap=argparse.ArgumentParser(description="NeuroThermo transition ensemble v1.1 reprojection")
    sp=ap.add_subparsers(dest="command",required=True)
    for cmd in ["validate","run"]:
        p=sp.add_parser(cmd);p.add_argument("--config",required=True)
    args=ap.parse_args();cfg=load_config(args.config)
    out=validate(cfg,args.config) if args.command=="validate" else run(cfg,args.config)
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=="__main__": main()
