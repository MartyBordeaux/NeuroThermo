# NeuroThermo dynamic characterization v2.1

This is the experimental-support-restricted successor to v2.0. It does **not** refit the Hindmarsh-Rose model. Frozen v3.9 cell parameters, accepted spiking sweeps, selected spikes, rheobase brackets, near-optimal identifiability alternatives, and recovered animal IDs are the inputs.

## Scientific change relative to v2.0

v2.0 standardized cells at fixed `lambda = J / J_rheo`. That was useful as a stress test, but many lambda values lay outside the experimentally sampled current range of individual cells. v2.1 removes that extrapolation from the primary analysis.

1. Model and experiment are evaluated only at the **actually recorded accepted spiking currents** of each primary cell.
2. Every near-optimal parameter alternative is compared with the best fit at the **same physical current J**. Alternative-specific rheobase is still recalculated, but it does not change the current used for the suprathreshold robustness comparison.
3. Cross-cell suprathreshold position is represented by

   `q = (J - J_rheo_best) / (J_max,observed - J_rheo_best)`.

   The default targets are `q = 0.25, 0.50, 0.75`. Interpolation is allowed only between actually observed spiking currents. Extrapolation is forbidden.
4. Scalar model metrics are restricted to the original v3.9 sweep support window (`0 .. fit_end_ms`) after exact first-spike alignment. Extra model spikes outside that window are not treated as supported phenotype.
5. Phase-resolved model profiles use only complete model cycles whose **aligned spike boundaries lie between the first and last selected experimental spikes**. The model state trajectory itself is never shifted or time-rescaled.
6. Experimental spike-train descriptors are interpolated to the same supported q targets, enabling model-versus-experiment comparison without synthetic current extrapolation.
7. Rheobase remains a separate absolute-excitability phenotype and is refined for the best fit and all near-optimal alternatives.

No genotype-level animal p-values are generated because the recovered animal count is insufficient for reliable population-level inference.

## Frozen cohort

Primary: 18 multi-sweep cells (12 WT, 6 SCA3). The two accepted single-sweep cells remain secondary and are excluded from the primary q analysis. `SCA3_05` is handled as an ordinary accepted primary cell.

## Run

```bash
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -e .
python3 -m dynamic_v2.cli validate --config configs/server_dynamic_v2_1.yaml
./run_server.sh configs/server_dynamic_v2_1.yaml
```

Outputs are described in `OUTPUTS.md`.
