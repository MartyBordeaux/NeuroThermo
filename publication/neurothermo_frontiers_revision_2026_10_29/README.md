# NeuroThermo Frontiers revision - 2026-10-29

## Status and scope

This is a complete working revision of the existing manuscript, not a new experimental study and not a journal submission. The supported result is input-protocol-dependent ordering of active-cycle HR state affinity and log-ISI affinity. The negative fixed-current result, failed original WT-WT criterion, limited endpoint-matched control overlap, and negative historical matched-noise/proxy checks remain explicit.

The primary observable excludes silent intervals and uses equal mass per accepted ISI cycle. It is not the stationary full-time occupancy used in the original manuscript. No new HR fitting was performed.

## Documents

- `main.tex` / `main.pdf`: revised manuscript with four main figures.
- `supplementary.tex` / `supplementary.pdf`: exact methods, provenance limits, historical evidence, two supplementary figures.
- `response_to_reviewers.tex` / `response_to_reviewers.pdf`: all 23 reviewer points, original comments, responses, and resolved section references.
- `references.bib`: separate bibliography; original database retained and four relevant primary references added.
- `AUTHOR_CHECKS.md`: unresolved archival information, not invented methods.

## Rebuild from the complete companion bundle

```bash
python code/verify_release.py
python code/render_figures.py
bash build.sh
```

The rendering script requires NumPy, pandas, SciPy, Numba, Matplotlib, and Pillow. A LaTeX installation with pdfLaTeX, BibTeX, natbib, and the supplied Frontiers class is needed. Build order matters: the response uses cross-references from the compiled manuscript and supplement. The original scientific analysis ZIPs target Python 3.9; the environment used here for figure rendering is separately recorded in `provenance/RENDER_ENVIRONMENT.json`.

## Data and code

`source_data/` contains authoritative server summaries for v1.1.0-v1.4.0, archival experimental endpoint values, frozen model inputs, native-audit images/provenance, and illustrative figure-source values. The three Fig. 2 trajectories replay the stored best-fit example; no new cohort result is inferred from them.

`analysis_packages/` in the complete companion bundle preserves the original four scientific packages. `provenance/SOURCE_ARCHIVES.json` identifies their hashes and those of the authoritative result archives. Full production simulation checkpoints remain in the named original result archives; compact tables sufficient for the manuscript summaries are included here. No raw experimental data or historical fit is overwritten.

## Publication branch

The revision is added to `publication/plos-v1.3`, retaining the historical release under `publication/neurothermo_plos_v1_2/`. The new directory is `publication/neurothermo_frontiers_revision_2026_10_29/`.

The connector-written Git revision records manuscript/supplement sources, bibliography, table/figure definitions, compact numerical summaries, and the numerical verification and rendering sources. Compiled PDFs, raster panels, complete source tables, original analysis ZIPs, and the complete point-by-point response are supplied in the downloadable companion bundle. This text-layer commit must not be described as an independent clean-clone reproduction of all new binary assets. Copy the complete companion bundle into this directory to obtain the exact full build environment.

## Interpretation and checks

`NUMERICAL_RECHECK.json` records independent recomputation of the manuscript numerators and medians. Technical verification does not imply a passed genotype-specificity test. Cell pairs, parameter combinations, routes, and repeated q levels are not independent animals. All complete-case and attempted denominators remain distinct.
