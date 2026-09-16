# v1.3.1 figure-only correction

The numerical v1.3 decomposition results are unchanged.

The three original `01_boundary_*.png` panels in v1.3.0 were blank because the plotting function filtered the pandas table with `b.mode`. In pandas, `DataFrame.mode` is a method, so that expression did not select the `mode` column and the plot received no rows.

v1.3.1 uses explicit column access (`b['mode']` and `b['stage']`) and regenerates:

- `01_boundary_WT_exit.png`
- `01_boundary_balance.png`
- `01_boundary_SCA3_entry.png`

No HR simulations were rerun and no numerical CSV values were changed.
