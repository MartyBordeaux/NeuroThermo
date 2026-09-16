from pathlib import Path
import argparse, json
from .core import run, validate

def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='cmd',required=True)
    for name in ['validate','run']:
        q=sub.add_parser(name)
        q.add_argument('--config',required=True)
    a=p.parse_args()
    cfg=Path(a.config)
    out=validate(cfg) if a.cmd=='validate' else run(cfg)
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
