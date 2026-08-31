#!/usr/bin/env python3
"""Reconstruct and verify the frozen calibration inputs from the publication bundle."""
from __future__ import annotations
import base64
import hashlib
import json
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
B64 = ROOT / "data" / "calibration_bundle" / "neurothermo_publication_calibration_inputs_v1.tar.gz.b64"
CAL = ROOT / "data" / "calibration"
ARCHIVE_SHA256 = "d08702f704da9930f5351f377eadcb9527b759f4de1ef31aacddc00f4b6d2b53"
EXPECTED = {
    "candidate_events_with_predictions.csv": "af35c327b313482f534aa59669a47e52a4078f912a5e342efcfddf0158455640",
    "frozen_accepted_spiking_sweeps_v3_5.csv": "dad46b831eb4613af4a49673f83854e4ef48b81d0934c087234562d81a447a54",
    "frozen_peak_overrides_v3_5.csv": "64e35808199e6108355b015b4ca9ded6070deed852927877e705ccf118e95069",
    "frozen_threshold_brackets_v3_5.csv": "47ba271e6b8d70704de1c49aaac3677c6ee21e3001f33faaafad8761177f9741",
    "frozen_v3_1_cell_fit_summary.csv": "85b1fa2c457e4affc0db438cf885b4406f61b943cbc08073fbdeb7f4b57f42f9",
    "frozen_v3_1_sweep_fit_summary.csv": "5663d59c35aeb105ee45b0c4c8606375210294f377a6ee3adcd771356a70ab12",
    "frozen_v3_1_identifiability.csv": "16e810e3331a0f6eb6bc1c815bb0e0d5574ee93966b95c00346522d5470957d1",
    "seed_cell_summary_v3_9.csv": "cb74bc0783c9fd1db11cacba13ccabd273cfc225e6dc019ab6e4215433dceb72",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def safe_extract(tf: tarfile.TarFile, target: Path) -> None:
    target = target.resolve()
    for member in tf.getmembers():
        out = (target / member.name).resolve()
        if target != out and target not in out.parents:
            raise RuntimeError(f"Unsafe archive member: {member.name}")
    tf.extractall(target)


def main() -> int:
    if not B64.exists():
        raise FileNotFoundError(B64)
    encoded = "".join(B64.read_text(encoding="ascii").split())
    raw = base64.b64decode(encoded, validate=True)
    got_archive = sha256_bytes(raw)
    if got_archive != ARCHIVE_SHA256:
        raise RuntimeError(f"Calibration archive SHA-256 mismatch: {got_archive} != {ARCHIVE_SHA256}")

    with tempfile.NamedTemporaryFile(suffix=".tar.gz") as tmp:
        tmp.write(raw); tmp.flush()
        with tarfile.open(tmp.name, "r:gz") as tf:
            safe_extract(tf, ROOT)

    errors = []
    for name, expected in EXPECTED.items():
        path = CAL / name
        if not path.exists():
            errors.append(f"missing {name}")
            continue
        got = sha256_file(path)
        if got != expected:
            errors.append(f"{name}: {got} != {expected}")
    if errors:
        raise RuntimeError("Calibration verification failed:\n" + "\n".join(errors))

    print(f"PASS calibration archive sha256={ARCHIVE_SHA256}")
    for name in sorted(EXPECTED):
        print(f"PASS {name} {EXPECTED[name]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
