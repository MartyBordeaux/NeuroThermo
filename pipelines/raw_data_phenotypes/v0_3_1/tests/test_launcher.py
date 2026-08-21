from __future__ import annotations

import os
from pathlib import Path
import subprocess


def test_launcher_ignores_stale_installed_package(tmp_path: Path) -> None:
    project_root = Path(__file__).resolve().parents[1]
    stale_root = tmp_path / "stale"
    stale_package = stale_root / "neurothermo_phenotypes"
    stale_package.mkdir(parents=True)
    (stale_package / "__init__.py").write_text('__version__ = "0.1.2"\n', encoding="utf-8")
    (stale_package / "__main__.py").write_text(
        'raise SystemExit("stale package was executed")\n', encoding="utf-8"
    )

    env = os.environ.copy()
    env["PYTHONPATH"] = str(stale_root)
    env["NEUROTHERMO_DATA_ROOT"] = str(tmp_path / "missing_data")
    completed = subprocess.run(
        ["bash", str(project_root / "run_neurothermo.sh")],
        cwd=project_root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 1
    assert "Using neurothermo-phenotypes 0.3.1" in completed.stdout
    assert str(project_root / "src" / "neurothermo_phenotypes") in completed.stdout
    assert "stale package was executed" not in completed.stderr
    assert "Required workbook not found" in completed.stderr
