# Validation contract and reference run

v2.1 preserves the frozen v3.9 primary cohort and enforces the following invariants:

- 18 primary multi-sweep cells: 12 WT and 6 SCA3.
- 113 accepted spiking sweeps and 4884 selected spikes in the complete accepted frozen set.
- 111 spiking sweeps and 4853 selected spikes in the primary multi-sweep set.
- q is defined from the **best-fit** refined rheobase and the highest actually observed accepted spiking current of the same cell.
- No q target is evaluated outside a cell's observed q range.
- Near-optimal alternatives are simulated at exactly the same physical J as the corresponding experimental sweep.
- Exact first-spike alignment is additive only; there is no time rescaling and no last-spike anchoring.
- Scalar model phenotype excludes aligned spikes outside the original `fit_end_ms` support window.
- Phase profiles use only complete model cycles whose aligned boundaries lie within the selected experimental first-to-last-spike interval.
- Rheobase refinement remains separate from suprathreshold q dynamics.
- Formal animal-level genotype p-values are not produced.

## Automated checks

`pytest`: 4 passed.

The tests verify frozen cohort counts, additive-shift invariance of ISI/train metrics, exact first-spike alignment plus support clipping, and strict prohibition of scalar extrapolation.

## Full reference run

The bundled frozen inputs completed a full run successfully.

- best-fit rheobase bracket pass: 18/18 cells;
- supported cells at q=0.25: 12/18;
- supported cells at q=0.50: 17/18;
- supported cells at q=0.75: 18/18;
- maximum number of distinct physical J values for the same cell/sweep across best and alternative solutions: 1;
- unsupported q rows with finite core phenotype metrics: 0;
- best-fit cells with phase profiles after q interpolation: 8 at q=0.25, 16 at q=0.50, 17 at q=0.75.

Diagnostic robustness in the reference run (20% threshold, not an inferential significance criterion):

- all evaluated same-current scalar combinations stable: 30.7%;
- all evaluated q-interpolated scalar combinations stable: 32.1%;
- same-current firing-rate combinations stable: 60.9%;
- q-interpolated firing-rate combinations stable: 62.9%;
- same-current mean-ISI combinations stable: 60.3%;
- q-interpolated mean-ISI combinations stable: 58.6%.

These values are reference diagnostics only and should be interpreted after inspecting the user-side run.
