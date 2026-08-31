# NeuroThermo Stage 1: QC-calibrated spike phenotype

This is a separate post-visual-QC branch. It does not modify the preceding
fixed-height analysis.

The detector is trained only on the event classes explicitly accepted during
visual QC. Absolute peak voltage is recorded for audit but deliberately omitted
from the classifier, so low-amplitude spikes such as those in SCA3_02 can be
recognized from morphology. Predictors are prominence, local-baseline amplitude,
rise slope, fall slope, half-width, and after-hyperpolarization.

The main cohort excludes WT_01 (bursting sensitivity cell), WT_02
(artifact-prone), WT_05 and WT_07 (no usable signal). SCA3_02 is included but is
never used to train the detector. Technically valid low-current sweeps remain
zero responses. Plateau-associated cessation is retained as a biological
outcome. Currents above 600 pA are outside the primary comparison.

Run:

```bash
source /root/venv/bin/activate
python3 -m pytest -q

nohup env \
  WT_ROOT=/root/neurothermo/WT \
  SCA3_ROOT=/root/neurothermo/SCA3 \
  OUTPUT=results_stage1_qc_calibrated \
  bash run_qc_calibrated.sh > stage1_qc_calibrated.log 2>&1 &
```

The primary endpoint is sustained firing over the full fixed stimulus interval.
An isolated spike remains in `total_spike_count` but not in
`sustained_spike_count`. Statistical outputs are repeated for the main cohort
and prespecified sensitivity cohorts.
