from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def make_trace(group: str, cell_index: int, current: float, time: np.ndarray, rng) -> np.ndarray:
    voltage = -68.0 + rng.normal(0.0, 0.35, len(time))
    active = (time >= 0.1) & (time <= 1.1)
    voltage[active] += 0.02 * current
    threshold = 100 if group == "WT" else 150
    rate = max(0.0, (current - threshold) * (0.055 if group == "WT" else 0.038))
    rate *= 1.0 + rng.normal(0, 0.05)
    if rate > 1:
        first = 0.14 + (0.012 if group == "SCA3" else 0.0)
        spike_times = np.arange(first, 1.1, 1.0 / rate)
        for spike in spike_times:
            width = 0.0008 if group == "WT" else 0.0011
            voltage += 95 * np.exp(-0.5 * ((time - spike) / width) ** 2)
            voltage -= 8 * np.exp(-0.5 * ((time - spike - 0.0025) / 0.0015) ** 2)
    return voltage


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(17)
    time = np.arange(0.0, 1.3, 0.0005)
    currents = [0, 100, 200, 300, 400, 500, 600]
    manifest = []
    for group in ["WT", "SCA3"]:
        for cell_index in range(4):
            cell_id = f"{group}_{cell_index + 1:02d}"
            rows = []
            capacitance = (400 + 20 * cell_index) if group == "WT" else (100 + 8 * cell_index)
            for sweep_index, current in enumerate(currents):
                command = np.where((time >= .1) & (time <= 1.1), current, 0.0)
                voltage = make_trace(group, cell_index, current, time, rng)
                rows.append(pd.DataFrame({"time_s": time, "voltage_mV": voltage, "current_trace_pA": command, "current_pA": current, "sweep_index": sweep_index}))
            path = out / f"{cell_id}.csv"
            pd.concat(rows, ignore_index=True).to_csv(path, index=False)
            manifest.append({"group": group, "cell_id": cell_id, "record_id": cell_id + "_CC", "path": path.name, "capacitance_pF": capacitance, "capacitance_10ms_pF": capacitance * .8, "capacitance_20ms_pF": capacitance, "capacitance_50ms_pF": capacitance * 1.2, "include": True})
    pd.DataFrame(manifest).to_csv(out / "manifest.csv", index=False)


if __name__ == "__main__":
    main()

