#!/usr/bin/env python3
"""Extract, verify, and provenance-record the exact frozen calibration inputs."""
from __future__ import annotations
import csv
import hashlib
import shutil
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data" / "calibration_bundle" / "neurothermo_publication_calibration_bundle_2026-08-31.tar.gz"
CAL = ROOT / "data" / "calibration"
ARCHIVE_SHA256 = "0c930506021826aec8ee2987fe83cd4a1537fa42b6d3fad335a5520fcbb610bd"

SOURCES = {
    "candidate_events_with_predictions.csv": "neurothermo/results_stage1_qc_fixed/00_calibration/candidate_events_with_predictions.csv",
    "frozen_accepted_spiking_sweeps_v3_5.csv": "neurothermo/v4/hr_cell_fit_v3_5/calibration/frozen_accepted_spiking_sweeps_v3_5.csv",
    "frozen_peak_overrides_v3_5.csv": "neurothermo/v4/hr_cell_fit_v3_5/calibration/frozen_peak_overrides_v3_5.csv",
    "frozen_threshold_brackets_v3_5.csv": "neurothermo/v4/hr_cell_fit_v3_5/calibration/frozen_threshold_brackets_v3_5.csv",
    "frozen_v3_1_cell_fit_summary.csv": "neurothermo/v4/hr_cell_fit_v3_9/calibration/frozen_v3_1_cell_fit_summary.csv",
    "frozen_v3_1_sweep_fit_summary.csv": "neurothermo/v4/hr_cell_fit_v3_9/calibration/frozen_v3_1_sweep_fit_summary.csv",
    "frozen_v3_1_identifiability.csv": "neurothermo/v4/hr_cell_fit_v3_9/calibration/frozen_v3_1_identifiability.csv",
    "seed_cell_summary_v3_9.csv": "neurothermo/v4/hr_cell_fit_v3_9/calibration/seed_cell_summary_v3_9.csv",
}
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
ROLES = {
    "candidate_events_with_predictions.csv": "automatic spike-candidate/QC classifier output upstream of frozen manual selection",
    "frozen_accepted_spiking_sweeps_v3_5.csv": "canonical accepted spiking-sweep manifest used by final HR fit",
    "frozen_peak_overrides_v3_5.csv": "canonical manual peak overrides used by final HR fit",
    "frozen_threshold_brackets_v3_5.csv": "canonical non-spiking/spiking rheobase brackets used as binary threshold constraint",
    "frozen_v3_1_cell_fit_summary.csv": "frozen pre-wide-refit cell-level baseline/provenance",
    "frozen_v3_1_sweep_fit_summary.csv": "frozen pre-wide-refit sweep-level baseline/provenance",
    "frozen_v3_1_identifiability.csv": "frozen pre-wide-refit identifiability baseline/provenance",
    "seed_cell_summary_v3_9.csv": "deterministic starting/seed cell summary for final v3.9 refit",
}
V36_EQUIVALENTS = {
    "frozen_accepted_spiking_sweeps_v3_5.csv": "neurothermo/v4/hr_cell_fit_v3_6/calibration/frozen_accepted_spiking_sweeps_v3_6.csv",
    "frozen_peak_overrides_v3_5.csv": "neurothermo/v4/hr_cell_fit_v3_6/calibration/frozen_peak_overrides_v3_6.csv",
    "frozen_threshold_brackets_v3_5.csv": "neurothermo/v4/hr_cell_fit_v3_6/calibration/frozen_threshold_brackets_v3_6.csv",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_stream(fobj) -> str:
    h = hashlib.sha256()
    for chunk in iter(lambda: fobj.read(1024 * 1024), b""):
        h.update(chunk)
    return h.hexdigest()


def csv_shape(path: Path) -> tuple[int, int]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return 0, 0
        return sum(1 for _ in reader), len(header)


def write_provenance() -> None:
    out = CAL / "CALIBRATION_PROVENANCE.tsv"
    fields = ["publication_name", "original_tar_member", "source_server_archive", "sha256", "rows", "columns", "destination", "role"]
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for name in SOURCES:
            rows, cols = csv_shape(CAL / name)
            w.writerow({
                "publication_name": name,
                "original_tar_member": SOURCES[name],
                "source_server_archive": "neurothermo_publication_calibration_bundle_2026-08-31.tar.gz",
                "sha256": EXPECTED[name],
                "rows": rows,
                "columns": cols,
                "destination": f"data/calibration/{name}",
                "role": ROLES[name],
            })


def main() -> int:
    if not ARCHIVE.exists():
        raise FileNotFoundError(ARCHIVE)
    got_archive = sha256_file(ARCHIVE)
    if got_archive != ARCHIVE_SHA256:
        raise RuntimeError(f"Calibration archive SHA-256 mismatch: {got_archive} != {ARCHIVE_SHA256}")
    CAL.mkdir(parents=True, exist_ok=True)
    with tarfile.open(ARCHIVE, "r:gz") as tf:
        names = set(tf.getnames())
        missing = [m for m in list(SOURCES.values()) + list(V36_EQUIVALENTS.values()) if m not in names]
        if missing:
            raise RuntimeError("Missing expected archive members:\n" + "\n".join(missing))
        for out_name, member_name in SOURCES.items():
            src = tf.extractfile(member_name)
            if src is None:
                raise RuntimeError(f"Could not read archive member: {member_name}")
            with (CAL / out_name).open("wb") as dst:
                shutil.copyfileobj(src, dst)
        for canonical, member_name in V36_EQUIVALENTS.items():
            src = tf.extractfile(member_name)
            if src is None:
                raise RuntimeError(f"Could not read archive member: {member_name}")
            v36_hash = sha256_stream(src)
            if v36_hash != EXPECTED[canonical]:
                raise RuntimeError(f"v3.5/v3.6 content mismatch for {canonical}: {v36_hash} != {EXPECTED[canonical]}")
    errors=[]
    for name, expected in EXPECTED.items():
        got = sha256_file(CAL/name) if (CAL/name).exists() else "MISSING"
        if got != expected:
            errors.append(f"{name}: {got} != {expected}")
    if errors:
        raise RuntimeError("Calibration verification failed:\n" + "\n".join(errors))
    write_provenance()
    print(f"PASS calibration archive sha256={ARCHIVE_SHA256}")
    for name in sorted(EXPECTED):
        print(f"PASS {name} {EXPECTED[name]}")
    print("PASS v3.5/v3.6 frozen selection files are byte-identical")
    print("PASS data/calibration/CALIBRATION_PROVENANCE.tsv")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
