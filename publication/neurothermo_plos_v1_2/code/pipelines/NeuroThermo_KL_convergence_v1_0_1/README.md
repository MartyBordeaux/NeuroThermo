# NeuroThermo KL convergence v1.0.1

This package decides whether the endpoint-relative KL result is sufficiently
stable to remain a principal result of the WT--SCA3 manuscript. It does not
assume that every support scenario has one universal crossing.

## Full design

- frozen core-secure cohort: 264 support scenarios in 32 dependent cell pairs;
- path grid: 31 coupled-path positions;
- integration steps: 0.05, 0.025, and 0.0125 ms;
- five frozen stochastic seeds;
- distributional views: full `xyz`, fast `xy`, and slow `z`;
- 122,760 pilot and 122,760 analysis simulations (245,520 total);
- a common scenario-specific density grid spanning the observed extrema of all
  15 step--seed paths;
- nested Brownian increments across integration steps;
- marker-first and curve-first aggregation;
- pair-, endpoint-cell-, animal-pair-, and leave-one-animal-out summaries.

The calculation is explicitly two-pass. The pilot pass stores only extrema,
then the analysis pass deterministically repeats the paths on the resulting
full-coverage grids. The v1.0.0 analysis pass took 1.13 h on the reference
2-CPU server, so v1.0.1 should take approximately 2.3 h. Both passes are
resumable after every complete scenario--step--seed path.

## Installation

```bash
unzip NeuroThermo_KL_convergence_v1_0_1.zip
cd NeuroThermo_KL_convergence_v1_0_1
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
./smoke_test.sh
```

## Full server run

The publication release contains the analysis-ready inputs at
`../../../data/inputs`, so no machine-specific input path is required:

```bash
./start_nohup.sh
tail -f kl_convergence_v1_0_1.log
```

Resume uses the identical command. A checkpoint is accepted only when its
scientific fingerprint matches the configuration, frozen input hashes, animal
mapping, and package version.

After completion:

```bash
./pack_results.sh
```

Return `neurothermo_kl_convergence_results_v1_0_1.zip`.

## Decision

Read `results_kl_convergence_v1_0_1/KL_CONVERGENCE_VERDICT.json`.

- `KEEP_AS_MAIN_RESULT`: all fatal and supporting gates passed.
- `KEEP_AS_ENSEMBLE_RESULT_WITH_LIMITATIONS`: all fatal gates for the primary
  full-state ensemble result passed, but at least one supporting sensitivity
  gate failed.
- `REMOVE_KL_RESULT`: at least one fatal primary-result gate failed.

A keep decision supports only an ensemble ordering along the constructed path.
It does not establish a universal scenario crossing, biological disease time,
causal precedence, irreversibility, or an early clinical biomarker.
