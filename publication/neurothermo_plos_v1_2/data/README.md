# Data inventory

All capacitances, current densities, and model parameters in this release are
stored in the final units used by the manuscript.

## Cohort provenance

`animal_to_cell_mapping.csv` is the authoritative cell-to-animal map.

- WT raw archive: 16 cells from 7 animal-days.
- WT fitted cohort: 13 cells from 6 animal-days.
- WT primary multi-sweep cohort: 12 cells from 6 animal-days.
- SCA3 raw archive: 9 cells from 2 animal-days (`DD20`, `DD24`).
- SCA3 fitted cohort: 7 cells from those same 2 animal-days.
- SCA3 primary multi-sweep cohort: 6 cells from those same 2 animal-days.

The recovered current-clamp archive is a subset of the broader experimental
material reported by Konno and colleagues; it is not treated as the complete
animal pool from that publication. The two SCA3 animal-days limit biological
replication, so the article does not report genotype-population inference.

## Data levels

- `inputs/` contains the 264 core-secure support scenarios and 32 dependent
  WT--SCA3 cell-pair stage definitions consumed by the heavy pipelines.
- `figure_source/` contains one tidy CSV per article panel or supporting table.
- `kl_convergence_v1_0_1/` contains complete aggregate full-coverage outputs.
  `scenario_markers.csv.gz` is compressed losslessly and is readable directly
  by pandas or R/readr.
- `nonequilibrium_geometry_v1_0_1/` contains the animal- and animal-pair-balanced
  outputs used for the nonequilibrium result. The animal-pair table is stored as
  lossless `csv.gz`. Large diagnostic-only estimator grids that are not used in
  the article are intentionally excluded.

Raw ABF recordings are not included in this repository release.
