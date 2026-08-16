# Cell-fit methodology — frozen v3.9

## Model

For each accepted cell, one Hindmarsh–Rose parameter vector `(b, r, s, kappa_I)` is fitted jointly across all accepted spiking current steps. The fixed model constants are `a=1`, `c=1`, `d=5`, and `x_R=-1.6`.

The final common parameter space is:

- `b: 0.5–7`
- `r: 1e-4–0.1`
- `s: 0.05–15`
- `kappa_I: 2e-4–2`

The v3.9 boundary stress changed only the lower `s` bound from 0.25 to 0.05. No cell reached the new lower boundary.

## Spike-train objective

Experimental spike times are the manually reviewed Stage-1 events. The cell objective uses Victor–Purpura spike-train distance plus a spike-count penalty. There is no separate ISI, adaptation, latency, or voltage-plateau loss.

For each spiking sweep the model train is shifted additively so that its first spike coincides with the experimental first spike. This removes absolute onset latency as a nuisance component while preserving every model ISI and the complete model spike count. Model spikes are not discarded after the shift.

The final spike is not anchored. No affine time rescaling is performed, because forcing both first and last spikes to coincide would remove information on train duration and average ISI.

## Rheobase constraint

For every cell, the nearest experimental silent/spiking bracket is retained. The model must satisfy `N_model(I0)=0` for the highest silent current and `N_model(I1)>=1` for the first spiking current. Non-spiking plateau voltage is not fitted.

## Identifiability

Practical identifiability is tested by searching for alternative parameter vectors that remain within the accepted loss tolerance while being separated from the optimum. Alternative solutions re-profile first-spike alignment. The separation criterion is not expanded in proportion to the wide v3.9 bounds, avoiding an artificial increase in PASS rate.

Final primary-cohort identifiability: all four parameters 5/18; `b` 11/18; `r` 6/18; `s` 9/18; `kappa_I` 10/18.
