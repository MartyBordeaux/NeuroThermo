# NeuroThermo Stage 1: fixed manual QC spike phenotype

This is a separate post-visual-QC branch. It does not modify the preceding
fixed-height analysis.

The morphology detector supplies the initial event decisions. The frozen
per-sweep decisions in `qc2.csv` are then applied as explicit overrides. Both
the algorithmic and final decisions are retained in the event table.

The primary `fixed_qc_all` cohort contains all 16 WT and 9 SCA3 cells reviewed
sweep by sweep. The earlier exclusions WT_01, WT_02, WT_05 and WT_07 define the
`conservative` sensitivity cohort. A second sensitivity cohort also excludes
SCA3_02. Plateau-associated cessation and valid zero responses are retained.
Currents above 600 pA are outside the primary comparison.

Run:

```bash
source /root/venv/bin/activate
python3 -m pytest -q

nohup env \
  WT_ROOT=/root/neurothermo/WT \
  SCA3_ROOT=/root/neurothermo/SCA3 \
  OUTPUT=results_stage1_qc_fixed \
  bash run_qc_fixed.sh > stage1_qc_fixed.log 2>&1 &
```

The primary endpoint is sustained firing over the full fixed stimulus interval.
An isolated spike remains in `total_spike_count` but not in
`sustained_spike_count`. Statistical outputs are repeated for the main cohort
and prespecified sensitivity cohorts.
