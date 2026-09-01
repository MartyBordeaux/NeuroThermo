#!/usr/bin/env python3
"""Generate or verify the SHA-256 manifest for the complete publication release tree."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "MANIFEST.sha256"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def entries() -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for p in sorted(x for x in ROOT.rglob("*") if x.is_file()):
        if p == MANIFEST:
            continue
        rel = p.relative_to(ROOT).as_posix()
        out.append((sha256_file(p), rel))
    return out


def render(rows: list[tuple[str, str]]) -> str:
    return "".join(f"{digest}  {rel}\n" for digest, rel in rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = ap.parse_args()

    current = render(entries())
    if args.write:
        MANIFEST.write_text(current, encoding="utf-8")
        print(f"WROTE {MANIFEST.relative_to(ROOT)} entries={current.count(chr(10))}")
        return 0

    if not MANIFEST.is_file():
        raise SystemExit("MANIFEST.sha256 is missing")
    frozen = MANIFEST.read_text(encoding="utf-8")
    if frozen != current:
        frozen_lines = set(frozen.splitlines())
        current_lines = set(current.splitlines())
        missing = sorted(frozen_lines - current_lines)[:20]
        new_or_changed = sorted(current_lines - frozen_lines)[:20]
        print("MANIFEST.sha256 mismatch")
        if missing:
            print("Missing/changed frozen entries:")
            print("\n".join(missing))
        if new_or_changed:
            print("New/changed current entries:")
            print("\n".join(new_or_changed))
        return 1
    print(f"MANIFEST_SHA256_PASS entries={current.count(chr(10))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
