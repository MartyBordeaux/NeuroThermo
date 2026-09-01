# Reproducibility and execution order

This release is organized around `code/run_full_analyses.sh`. All publication commands below are intended to be run from a clean clone and use repository-relative paths. The portable Python 3.9.25 publication environment is frozen in `environment/requirements-python39-historical.txt` and is referenced by the release-level `requirements.txt`. The filename is retained for provenance compatibility.

## 1. Environment

From `publication/neurothermo_plos_v1_2/`:

```bash
python3.9 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -e code/pipelines/NeuroThermo_endpoint_ensemble_v1_0_1 --no-deps
```

The publication-tested core versions are Python 3.9.25, NumPy 1.25.2, pandas 2.2.3, SciPy 1.13.1, scikit-learn 1.6.1, pyABF 2.3.8, Matplotlib 3.9.4, Numba 0.60.0, llvmlite 0.43.0, and statsmodels 0.14.6. NumPy 1.25.2 is a deliberate portability pin: with identical decoded ABF arrays, NumPy 1.26.4 produced CPU-dependent candidate identities on hosted runners, while 1.25.2 reproduced all 6800 frozen candidate identities with zero missing or extra events across the tested SciPy 1.9.3--1.13.1 matrix. The lock file contains the complete portable environment.

## 2. Release integrity and frozen inputs

```bash
bash code/run_full_analyses.sh raw-integrity
bash code/run_full_analyses.sh prepare
bash code/run_full_analyses.sh prepare-upstream
bash code/run_full_analyses.sh transition-frozen
bash code/run_full_analyses.sh transition-integrity
python code/generate_release_manifest.py --check
```

`MANIFEST.sha256` is the whole-release file manifest. It contains every regular file under the release root except `MANIFEST.sha256` itself. The calibration preparation step verifies the immutable calibration archive, materializes the canonical frozen CSV inputs, and creates `data/calibration/CALIBRATION_PROVENANCE.tsv`.

The manifest must be regenerated only after all release content, code and documentation changes are complete. A successful earlier manifest check is not evidence for a later commit that changes the release tree.

## 3. Raw ABF to fixed spike/QC layer

The automated production front end is:

```bash
bash code/run_full_analyses.sh qc-tests
bash code/run_full_analyses.sh qc-recompute
```

The fixed-QC pipeline reads the 50 committed current-clamp ABF recordings, extracts candidate peaks, calibrates the classifier under the frozen rules, and applies the frozen visual-QC decisions in `qc2.csv`.

Clean-clone validation under the pinned historical environment reproduces the frozen event/QC layer exactly at the scientific-decision level:

- 6800 candidate events in both the frozen and recomputed tables;
- 6217 classifier-positive (`algorithm_detected`) events in both tables;
- 6039 final `fixed_qc_detected` events in both tables;
- 186 events changed by fixed visual QC in both tables;
- zero missing or extra candidate identities;
- zero decision mismatches.

The clean-clone gate enforces this layer rather than documenting a single favorable runner: it requires all 6800 candidate identities, 6217 classifier-positive events, 6039 final fixed-QC detections, 186 visual-QC changes, zero missing or extra identities, zero decision mismatches, and `spike_probability` agreement within `rtol=1e-12, atol=1e-12`. A release commit is valid only when this gate and the figure-source and transition-chain gates all pass on that same SHA.

The manual visual-audit program is retained in `code/pipelines/NeuroThermo_spike_visual_qc/`. The final audited selections themselves are immutable publication inputs (`frozen_accepted_spiking_sweeps_v3_5.csv`, `frozen_peak_overrides_v3_5.csv`, and `frozen_threshold_brackets_v3_5.csv`) and do not require a reviewer to repeat an interactive manual audit.

## 4. Restricted Hindmarsh-Rose fit and characterization

Validate the exact frozen v3.9 implementation and inputs:

```bash
bash code/run_full_analyses.sh cellfit-validate
```

Full four-parameter refit (`b`, `r`, `s`, `kappa_I`) and identifiability runs are exposed separately because they are computationally expensive:

```bash
bash code/run_full_analyses.sh cellfit
bash code/run_full_analyses.sh cellfit-identify
```

The fit uses spike times, exact first-spike additive alignment, the frozen accepted spiking sweeps/manual peak overrides, and the binary rheobase constraint. Non-spiking sweeps are used only for the threshold bracket. The last spike is not fixed and time is not rescaled between the first and last spike.

Post-fit characterization is reproduced with:

```bash
bash code/run_full_analyses.sh characterization
```

The compatibility XLSX required by the historical characterization script is generated deterministically from the canonical `data/animal_id_recovery/accepted_cohort.csv`; the spreadsheet is not a second scientific source.

## 5. Dynamic and endpoint ensemble layers

```bash
bash code/run_full_analyses.sh dynamic-validate
bash code/run_full_analyses.sh dynamic
bash code/run_full_analyses.sh endpoint-validate
bash code/run_full_analyses.sh endpoint
```

The committed exact/frozen upstream layers allow immediate validation and downstream reproduction. Recomputed outputs are written under `results/recomputed/` and do not overwrite frozen publication results.

The validated v3.9 layer contains 20 accepted cells, 18 primary multi-sweep cells, 113 spiking fit sweeps and 4884 selected spikes after overrides. The primary dynamic cohort contains 18 cells (12 WT, 6 SCA3), 111 spiking sweeps and 4853 selected spikes.

## 6. Transition chain

The executable transition DAG is v1.0 -> v1.1 -> v1.2 -> v1.2.1 -> v1.3.

```bash
bash code/run_full_analyses.sh transition-frozen
bash code/run_full_analyses.sh transition-v1-0-validate
bash code/run_full_analyses.sh transition-v1-1-validate
bash code/run_full_analyses.sh transition-v1-1
bash code/run_full_analyses.sh prepare-transition-v1-2
bash code/run_full_analyses.sh transition-v1-2-validate
bash code/run_full_analyses.sh transition-v1-2-1-validate-frozen
bash code/run_full_analyses.sh transition-v1-3-validate-frozen
bash code/run_full_analyses.sh transition-integrity
```

Full simulation commands are:

```bash
bash code/run_full_analyses.sh transition-v1-0
bash code/run_full_analyses.sh transition-v1-2
bash code/run_full_analyses.sh transition-v1-2-1
bash code/run_full_analyses.sh transition-v1-3
```

v1.1 is a zero-new-simulation reprojection stage. The v1.1 -> v1.2 assembly is explicit and verifies the recalculated staging inputs against the frozen references. v1.2 and v1.3 are the heavy HR-state stages and use committed checkpoint/result trees for routine clean-clone validation; their full commands support the historical checkpoint/resume workflow.

The constructed path coordinate `p` orders model states. It is not interpreted as disease time, a causal trajectory, or evidence of irreversible one-way progression.

## 7. KL, nonequilibrium geometry, Fisher check, and figures

Lightweight publication validation of the expensive stochastic analyses is:

```bash
bash code/run_full_analyses.sh kl-validate
bash code/run_full_analyses.sh nonequilibrium-validate
```

Full stochastic recomputation remains available as:

```bash
bash code/run_full_analyses.sh kl
bash code/run_full_analyses.sh nonequilibrium
```

The reviewer-facing deterministic figure replay is:

```bash
bash code/reproduce_figures.sh
```

This command rebuilds deterministic publication source tables, deletes all existing `Fig*.pdf` and `Fig*.png` renders, regenerates the canonical figures, and requires exactly the following eight stems in both PDF and PNG formats (16 files total):

- `Fig1_endpoint_phenotype`
- `Fig2_transition_staging`
- `Fig3_intrinsic_drive_decomposition`
- `Fig4_thermodynamic_information_geometry`
- `Fig5_nonequilibrium_geometry`
- `FigS1_support_restricted_dynamics`
- `FigS2_multiseed_marker_robustness`
- `FigS3_KL_full_coverage_convergence`

`Fig4_thermodynamic_information_geometry` has one canonical implementation: `render_figures.py` delegates to `render_fig4_multiseed.py`; the master entrypoint does not call the multiseed renderer a second time. The historical `figS2()` thermodynamic-robustness function retained inside `render_figures.py` is not invoked by the canonical entrypoint and is not part of the frozen 16-file release contract.

The Fisher-information-related local consistency output is part of the nonequilibrium-geometry analysis (`local_kl_fisher_check.csv`); it is not a separate standalone Fisher pipeline in the frozen source inventory.

The current manuscript result-to-code map does not use the historical pathwise-temporal-order package or the historical PI/Fourier package. They are therefore explicitly outside the publication-critical DAG for this release and are not required for reviewer clean-clone reproduction. They should be added only if a manuscript result is restored that explicitly depends on them. The manuscript's KL result is covered by `NeuroThermo_KL_convergence_v1_0_1`, including its full-state and marginal analyses; no separate publication-critical 2D-KL package is required for the current result map.

## 8. Repository footprint

The publication directory is intentionally large because it contains raw recordings plus exact frozen/checkpoint result layers needed to avoid forcing reviewers to rerun the computationally expensive HR-state and stochastic analyses. At the final pre-release audit it is approximately 945 MB (about 901 MiB).

The largest tracked file, `transition_surface_states.csv.gz`, is approximately 99.3 MB (about 94.7 MiB), close to GitHub's per-file limit. The clean-clone preflight therefore rejects any publication file larger than 95 MiB. The current file remains below that stricter release guard, but this margin should not be reduced by later edits. If the public snapshot cannot be pushed reliably with the current payload, the preferred remedy is a versioned split/reassembly or archival-data strategy rather than silently dropping the frozen state layer.

## 9. Acceptance criterion

A publication snapshot is ready to expose publicly only when, on the same final commit:

1. the whole-release `MANIFEST.sha256` verifies;
2. release hygiene contains no tracked runtime caches, `.pyc`, or historical `.log` files;
3. the main clean-clone preflight passes, including exact raw-to-QC scientific decisions;
4. the figure-source preflight passes and proves fresh generation of the exact 16-file figure contract;
5. the transition-chain preflight passes, including v1.1 reprojection and v1.1 -> v1.2 input assembly;
6. every manuscript result is mapped to a committed script, its inputs, and its frozen or recomputed outputs;
7. a final immutable public snapshot/tag/release is created from the validated commit.

Historical successful CI runs on earlier commits are provenance, not evidence that the current head passes these gates.
