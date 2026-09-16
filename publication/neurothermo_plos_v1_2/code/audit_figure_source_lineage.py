#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'data' / 'figure_source'
UPSTREAM_DIRS = [
    ROOT / 'data' / 'endpoint_ensemble_v1_0_results',
    ROOT / 'data' / 'transition_v1_1_results',
    ROOT / 'data' / 'transition_v1_2_1_results',
    ROOT / 'data' / 'transition_v1_3_results',
]
OUT = ROOT / 'docs' / 'FIGURE_SOURCE_LINEAGE_AUDIT.json'


def read_table(path: Path):
    try:
        if path.suffix == '.gz' or path.name.endswith('.csv.gz'):
            return pd.read_csv(path, compression='gzip')
        return pd.read_csv(path)
    except Exception:
        return None


def canon(df: pd.DataFrame) -> pd.DataFrame:
    z = df.copy()
    for c in z.columns:
        if pd.api.types.is_numeric_dtype(z[c]):
            z[c] = pd.to_numeric(z[c], errors='coerce').round(12)
        else:
            z[c] = z[c].astype(str)
    return z


def projection_match(source: pd.DataFrame, candidate: pd.DataFrame):
    if not set(source.columns).issubset(candidate.columns):
        return None
    s = canon(source[source.columns])
    c = canon(candidate[source.columns])
    # Compare row multisets by stable string keys; order independent.
    sk = s.astype(str).agg('\x1f'.join, axis=1).value_counts().sort_index()
    ck = c.astype(str).agg('\x1f'.join, axis=1).value_counts().sort_index()
    common = sk.index.intersection(ck.index)
    matched = int(sum(min(int(sk[k]), int(ck[k])) for k in common))
    return {
        'source_rows': int(len(s)),
        'candidate_rows': int(len(c)),
        'matched_source_rows': matched,
        'source_row_fraction': float(matched / len(s)) if len(s) else 1.0,
        'exact_projection_multiset': bool(matched == len(s)),
    }


def main():
    candidates = []
    for d in UPSTREAM_DIRS:
        if not d.is_dir():
            continue
        for p in sorted(list(d.glob('*.csv')) + list(d.glob('*.csv.gz'))):
            df = read_table(p)
            if df is None:
                continue
            candidates.append((p, df))

    report = {'figure_sources': {}, 'candidate_count': len(candidates)}
    for sp in sorted(FIG.glob('fig[123]_*.csv')):
        s = pd.read_csv(sp)
        entries = []
        for cp, c in candidates:
            shared = [x for x in s.columns if x in c.columns]
            rec = {
                'candidate': str(cp.relative_to(ROOT)),
                'candidate_rows': int(len(c)),
                'source_rows': int(len(s)),
                'source_columns': list(s.columns),
                'candidate_columns': list(c.columns),
                'shared_columns': shared,
                'shared_column_count': len(shared),
                'source_column_count': len(s.columns),
            }
            pm = projection_match(s, c)
            if pm is not None:
                rec['projection_match'] = pm
            if len(shared) >= max(2, int(0.5 * len(s.columns))) or pm is not None:
                entries.append(rec)
        entries.sort(key=lambda x: (
            -x.get('projection_match', {}).get('source_row_fraction', -1),
            -x['shared_column_count'], x['candidate']
        ))
        report['figure_sources'][sp.name] = {
            'rows': int(len(s)),
            'columns': list(s.columns),
            'top_candidates': entries[:12],
        }

    OUT.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
