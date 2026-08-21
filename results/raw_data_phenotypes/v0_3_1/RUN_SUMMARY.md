# NeuroThermo model-free phenotype run

## Dataset

| group | cells | analysed_sweeps | qc_pass_clean | qc_warning | qc_fatal | spiking | thermo_eligible | protocol_excluded |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SCA3 | 7 | 91 | 52 | 39 | 0 | 33 | 29 | 0 |
| WT | 13 | 169 | 130 | 39 | 0 | 67 | 63 | 14 |

## Primary interpretation

- Current in pA is the matched experimental protocol axis.
- J = I/Cm is emitted only inside shared WT/SCA3 support, including 10/20/50 ms capacitance sensitivity.
- Non-spiking sweeps enter rheobase brackets; thermodynamic trace metrics require at least the configured number of spikes.
- Curated Stage-1 events replace automatic peak detection in production.
- Production spike events are restricted to the hash-checked frozen v3.5 sweep keys; every rejected or window-excluded event is audited.
- Latency is preserved and never removed by trace alignment.
- Inference is cell-level. Sweep fragments are not treated as independent biological replicates.
- Integrated inference is restricted to unconditional features with complete configured common-current support.
- Conditional high-current inference is two-part: all accepted cells enter binary availability, while values are compared only where physically defined.
- No cell is removed globally for an incomplete conditional curve and no value is imputed on a non-spiking sweep.
- The disease coordinate gives equal weight to log-capacitance, mean-ISI and predictive-information domains after robust WT-reference scaling.
- q=0 and q=1 are the WT and SCA3 endpoint medians. q is not a disease probability and injected current is not disease time.
- Path-KL excess is the sole primary irreversibility metric; raw and signed bias-corrected variants are diagnostics.

## Scope limits

External work is clamp-supplied incremental electrical work. Path KL is a bias-controlled observable irreversibility proxy. Information entropies are not heat or metabolic entropy. This endpoint analysis separates WT and SCA3 phenotypes but does not identify biological disease time or a point of no return.

Current-resolved tests: 123.
Response-curve global tests: 15.
Cell-scalar tests: 12.
Integrated cell-level tests: 6.
All-cell two-part tests: 4.
Protocol-excluded sweeps: 14.

## Model-free disease coordinate

- Descriptive WT/SCA3 endpoint AUC: 1.0000.
- Internal leave-one-WT-out AUC: 1.0000; minimum SCA3-minus-WT margin: 0.0982384.
- Exact recomputed two-sided permutation p: 0.00230908.
- SCA3 cells above the robust WT boundary: 7/7.
- Both AUC values are internal diagnostics, not external or clinical performance.
- WT-exit is an operational multi-domain marker. A biological transition time requires intermediate or longitudinal samples.

## Curated event audit

- SCA3 / excluded_not_frozen_sweep: 257 events.
- WT / excluded_not_frozen_sweep: 890 events.
- WT / excluded_outside_common_current_grid: 774 events.
- SCA3 / excluded_peak_override: 5 events.
- WT / excluded_peak_override: 2 events.
- SCA3 / used_frozen_event: 680 events.
- WT / used_frozen_event: 3431 events.
