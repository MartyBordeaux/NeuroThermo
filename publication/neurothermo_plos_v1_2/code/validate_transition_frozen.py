#!/usr/bin/env python3
from pathlib import Path
import hashlib, sys

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
EXPECTED = {
    DATA/'transition_v1_0_frozen'/'cell_q75_protocol_anchors.csv': 'fe986b58f072ef0f00d252ca799e1f96b6e5b998a9acc3eade6a78a7967ab136',
    DATA/'transition_v1_2_frozen'/'cell_q75_protocol_anchors.csv': 'fe986b58f072ef0f00d252ca799e1f96b6e5b998a9acc3eade6a78a7967ab136',
    DATA/'transition_v1_3_frozen'/'cell_q75_protocol_anchors.csv': 'fe986b58f072ef0f00d252ca799e1f96b6e5b998a9acc3eade6a78a7967ab136',
    DATA/'transition_v1_1_frozen'/'q75_reference_cells.csv': 'd2b56731b8138fbc2106dcdd04d4b3d0cc467f9846e4d86cd9c21ea83b76a6f3',
}

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

for p, expected in EXPECTED.items():
    if not p.is_file():
        print(f'FAIL missing {p.relative_to(ROOT)}', file=sys.stderr)
        raise SystemExit(1)
    got = sha256(p)
    if got != expected:
        print(f'FAIL SHA mismatch {p.relative_to(ROOT)}: {got} != {expected}', file=sys.stderr)
        raise SystemExit(1)
    print(f'PASS {p.relative_to(ROOT)} {got}')

# The q75 anchors were frozen independently into v1.0, v1.2 and v1.3.
# Identity across all three copies is itself a provenance check.
anchors = [p for p in EXPECTED if p.name == 'cell_q75_protocol_anchors.csv']
raw = [p.read_bytes() for p in anchors]
if not (raw[0] == raw[1] == raw[2]):
    print('FAIL q75 anchor copies are not byte-identical', file=sys.stderr)
    raise SystemExit(1)
print('PASS v1.0/v1.2/v1.3 q75 anchor copies are byte-identical')
print('TRANSITION_FROZEN_INTEGRITY_PASS')
