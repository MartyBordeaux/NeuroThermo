# v1.0.1

- Made the coarse Markov time-reversal test the sole formalism decision basis.
- Marked the continuous-current branch diagnostic-invalid when divergence QC
  fails; its circulation and EPR-like values are never used as NESS evidence.
- Centered every path score before calculating Fisher information.
- Derived friction amplitude from the centered density-based FI and correlation
  time from the sampled force trajectory.
- Added Fisher chain-rule and force-variance QC fields.
- Replaced repeated nearest-neighbour adaptive schedules with strictly
  increasing 15-state schedules and matched the linear schedule to 15 states.
- Added protocol-performance and numerical-validation verdict files.
- Added scientific-fingerprint-protected compressed Markov caches.
- Added scenario/pair/seed provenance to all cycle affinities.
- Added endpoint-group membership and oscillatory-endpoint sensitivity outputs.
- Added the recovered cell-to-animal mapping, four-WT/two-SCA3 validation, and
  animal-pair-balanced geometry, protocol, and fluctuation summaries.
- Reworked diagnostic figures: no titles, logarithmic FI/friction axes, no
  invalid continuous-current evidence panel, and a visible log-path-ratio panel.
