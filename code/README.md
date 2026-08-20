# Readable working code

This subtree contains the working source used for frozen NeuroThermo analysis stages. Source is stored as normal, browsable Python/shell/YAML files rather than only as result tables or opaque archives.

The repository also contains results and documentation from parallel NeuroThermo work. Those existing paths are intentionally not modified by this code freeze.

Current readable snapshots in this freeze include the final HR cell-fit stage, post-fit characterization, the transition decomposition stage, and thermodynamic transition v1.0.1. Additional historical intermediate transition snapshots can be added under separate versioned directories without changing parallel-work paths.

Bulky numerical checkpoints such as `.npz` files are not part of the GitHub freeze; they are reproducibility caches rather than human-readable source or compact scientific results.
