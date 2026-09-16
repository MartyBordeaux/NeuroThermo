from __future__ import annotations
import argparse, json, yaml
from pathlib import Path
from .core import validate, run_all

def load_cfg(path):
    p=Path(path).resolve(); cfg=yaml.safe_load(p.read_text(encoding='utf-8'))
    base=p.parent.parent
    for key in ['directory','archive']:
        val=cfg.get('input',{}).get(key)
        if val and not Path(val).is_absolute(): cfg['input'][key]=str((base/val).resolve())
    val=cfg.get('output',{}).get('dir')
    if val and not Path(val).is_absolute(): cfg['output']['dir']=str((base/val).resolve())
    return cfg

def main():
    ap=argparse.ArgumentParser(); sub=ap.add_subparsers(dest='cmd',required=True)
    for name in ['validate','run']:
        p=sub.add_parser(name);p.add_argument('--config',required=True)
    args=ap.parse_args(); cfg=load_cfg(args.config)
    res=validate(cfg) if args.cmd=='validate' else run_all(cfg)
    print(json.dumps(res,indent=2,sort_keys=True))
if __name__=='__main__': main()
