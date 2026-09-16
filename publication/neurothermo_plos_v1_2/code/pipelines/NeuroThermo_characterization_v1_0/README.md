# NeuroThermo post-fit characterization v1.0.0

This stage characterizes the final v3.9 cell fits and merges recovered animal identifiers.

Primary analysis set: accepted multi-sweep cells. Single-sweep cells are retained as secondary descriptive records only.

The script deliberately does not report formal group p-values because the current recovered animal structure contains only two animals per group. Cell-level contrasts are descriptive effect sizes for the recorded-cell ensemble; animal-level summaries are medians per recovered animal.

Run:

```bash
python3 run_characterization.py \
  --results results_cellfit_v3_9.zip \
  --animal-map NeuroThermo_animal_id_recovery.xlsx \
  --outdir results_characterization_v1_0
```

Outputs include cell- and animal-level master tables, descriptive group effect sizes, identifiability summaries, Spearman association tables, publication-oriented figures, an Excel workbook, and REPORT.md.
