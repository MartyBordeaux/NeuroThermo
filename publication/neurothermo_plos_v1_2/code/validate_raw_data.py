#!/usr/bin/env python3
from pathlib import Path
import csv, hashlib, sys

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data' / 'raw'
MANIFEST = ROOT / 'docs' / 'RAW_DATA_MANIFEST.tsv'
PROV = ROOT / 'data' / 'upstream_server_bundle' / 'RAW_TRANSITION_FROZEN_ARCHIVE_SHA256.txt'
EXPECTED_ARCHIVE_SHA = '53e3384e9ce110a91fa971046a298e7cca2446ea7e8dd7185821ee95c3cbb43f'
EXPECTED_COUNTS = {'WT': 32, 'SCA3': 18}


def fail(msg):
    print('FAIL ' + msg, file=sys.stderr)
    raise SystemExit(1)


def ok(msg):
    print('PASS ' + msg)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

if not PROV.exists():
    fail('missing raw/transition archive provenance')
prov_sha = PROV.read_text(encoding='utf-8').split()[0]
if prov_sha != EXPECTED_ARCHIVE_SHA:
    fail(f'archive SHA mismatch: {prov_sha} != {EXPECTED_ARCHIVE_SHA}')
ok(f'raw/transition archive SHA256={prov_sha}')

if not MANIFEST.exists():
    fail('missing docs/RAW_DATA_MANIFEST.tsv')

with MANIFEST.open(encoding='utf-8', newline='') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))

if len(rows) != 50:
    fail(f'manifest row count {len(rows)} != 50')
ok('manifest rows=50')

seen = set()
counts = {'WT': 0, 'SCA3': 0}
total_bytes = 0
for r in rows:
    group = r['group']
    if group not in counts:
        fail(f'unknown group {group!r}')
    rel = r['relative_path']
    p = ROOT / rel
    key = (group, r['filename'])
    if key in seen:
        fail(f'duplicate manifest entry {key}')
    seen.add(key)
    if not p.is_file():
        fail(f'missing raw file {rel}')
    if p.suffix.lower() != '.abf':
        fail(f'non-ABF raw file in manifest: {rel}')
    if p.name != r['filename']:
        fail(f'filename/path mismatch for {rel}')
    size = p.stat().st_size
    if size != int(r['size_bytes']):
        fail(f'size mismatch {rel}: {size} != {r["size_bytes"]}')
    got = sha256(p)
    if got != r['sha256']:
        fail(f'SHA mismatch {rel}: {got} != {r["sha256"]}')
    counts[group] += 1
    total_bytes += size

for group, expected in EXPECTED_COUNTS.items():
    if counts[group] != expected:
        fail(f'{group} count {counts[group]} != {expected}')
    disk = sorted((RAW / group).glob('*.abf'))
    if len(disk) != expected:
        fail(f'{group} filesystem ABF count {len(disk)} != {expected}')
    ok(f'{group} ABF files={expected}')

all_disk = list(RAW.rglob('*.abf'))
if len(all_disk) != 50:
    fail(f'filesystem ABF count {len(all_disk)} != 50')

ok(f'total raw bytes={total_bytes}')
print('RAW_DATA_INTEGRITY_PASS')
