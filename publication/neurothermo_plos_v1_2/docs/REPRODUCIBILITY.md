# Reproducibility and execution order

This release is organized around `code/run_full_analyses.sh`. All publication commands below are intended to be run from a clean clone and use repository-relative paths. The historical Python 3.9.25 numerical environment is frozen in `environment/requirements-python39-historical.txt` and is referenced by the release-level `requirements.txt`.

## 1. Environment

From `publication/neurothermo_plos_v1_2/`:

```bash
python3.9 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -e code/pipelines/NeuroThermo_endpoint_ensemble_v1_0_1 --no-deps
```

The publication-tested core versions are Python 3.9.25, NumPy 1.26.4, pandas 2.2.3, SciPy 1.13.1, scikit-learn 1.6.1, pyABF 2.3.8, Matplotlib 3.9.4, Numba 0.60.0, llvmlite 0.43.0, and statsmodels 0.14.6. The lock file contains the complete captured environment.

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

## 3. Raw ABF to fixed spike/QC layer

The automated production front end is:

```bash
bash code/run_full_analyses.sh qc-tests
bash code/run_full_analyses.sh qc-recompute
```

The fixed-QC pipeline reads the 50 committed current-clamp ABF recordings, extracts candidate peaks, calibrates the frozen classifier rules, and applies the frozen visual-QC decisions in `qc2.csv`. The historical reference contains 6800 candidate events and 6039 final `fixed_qc_detected` events.

Current release status: raw-file integrity, fixed-QC unit tests, and independent candidate-identity audits pass. A hosted full-run consistency check has intermittently produced 6799 candidate rows / 6046 final detections instead of the frozen 6800 / 6039, although a dedicated threshold-boundary audit reproduced all 6800 frozen candidate identities with no missing or extra peaks and found no candidate within numerical epsilon of the frozen peak thresholds. This discrepancy is therefore treated as an unresolved execution/reproducibility issue. The publication gate must remain red until the full run reproduces the frozen candidate identities and final decisions; scientific thresholds are not relaxed to make the check pass.

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

The fit uses spike times, exact first-spike additive alignment, the frozen accepted spiking sweeps/manual peak overrides, and the binary rheobase constraint. Non-spiking sweeps are used only for the threshold bracket.

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

The committed frozen layers allow immediate validation and figure/source-table reproduction; recomputation writes under `results/recomputed/` and does not overwrite frozen publication results.

## 6. Transition chain

The executable transition DAG is v1.0 -> v1.1 -> v1.2 -> v1.2.1 -> v1.3.

```bash
bash code/run_full_analyses.sh transition-v1-0-validate
bash code/run_full_analyses.sh transition-v1-1-validate
bash code/run_full_analyses.sh transition-v1-1
bash code/run_full_analyses.sh prepare-transition-v1-2
bash code/run_full_analyses.sh transition-v1-2-validate
bash code/run_full_analyses.sh transition-v1-2-1-validate-frozen
bash code/run_full_analyses.sh transition-v1-3-validate-frozen
```

Full simulation commands are:

```bash
bash code/run_full_analyses.sh transition-v1-0
bash code/run_full_analyses.sh transition-v1-2
bash code/run_full_analyses.sh transition-v1-2-1
bash code/run_full_analyses.sh transition-v1-3
```

v1.1 is a zero-new-simulation reprojection stage. The v1.1 -> v1.2 assembly is explicit and verifies the recalculated staging inputs against the frozen references. v1.2 and v1.3 are the heavy HR-state stages and use committed checkpoint/result trees for routine clean-clone validation; their full commands support the historical checkpoint/resume workflow.

## 7. KL, nonequilibrium geometry, Fisher check, and figures

```bash
bash code/run_full_analyses.sh kl
bash code/run_full_analyses.sh nonequilibrium
bash code/run_full_analyses.sh figure-source
bash code/run_full_analyses.sh figures-python
```

The Fisher-information-related local consistency output is part of the nonequilibrium-geometry analysis (`local_kl_fisher_check.csv`); it is not a separate standalone Fisher pipeline in the frozen source inventory.

The pathwise temporal-order and frozen PI/Fourier repeat source trees are being incorporated as separate publication stages from their exact server snapshots. Until that import is complete, they must not be represented as clean-clone executable stages in this document or the master runner.

## 8. Acceptance criterion

A publication snapshot is ready to expose publicly only when, on the same final commit:

1. the whole-release `MANIFEST.sha256` verifies;
2. release hygiene contains no tracked runtime caches, `.pyc`, or historical `.log` files;
3. the main clean-clone preflight passes, including exact raw-to-QC consistency;
4. the figure-source preflight passes;
5. the transition-chain preflight passes;
6. every manuscript result is mapped to a committed script, its inputs, and its frozen or recomputed outputs.

Historical successful CI runs on earlier commits are provenance, not evidence that the current head passes these gates.
