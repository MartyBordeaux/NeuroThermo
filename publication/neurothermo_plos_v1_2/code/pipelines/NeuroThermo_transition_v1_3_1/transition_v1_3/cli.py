from __future__ import annotations
import argparse
from .decomposition import validate, run

def main():
    p=argparse.ArgumentParser(prog='transition_v1_3')
    sub=p.add_subparsers(dest='cmd',required=True)
    for name in ['validate','run']:
        sp=sub.add_parser(name); sp.add_argument('--config',required=True)
    a=p.parse_args()
    if a.cmd=='validate': validate(a.config)
    else: run(a.config)

if __name__=='__main__': main()
