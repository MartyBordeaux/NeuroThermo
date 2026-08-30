# NeuroThermo nonequilibrium geometry v1.0.1

This is a standalone, resumable server pipeline for testing whether the frozen
WT-to-SCA3 Hindmarsh–Rose path admits an equilibrium-potential interpretation.
It does **not** assume that Fisher information creates a Hamiltonian. Instead it
tests that assumption before allowing classical Jarzynski–Crooks language.

## What it calculates

1. Deterministic preflight at both endpoints: all fixed points, Jacobian
   eigenvalues, oscillation classification, and continuation in applied current.
2. Stationary density and nonequilibrium potential
   `phi(x,y,z;p) = -log rho_ss(x,y,z;p)` on one common grid per scenario/seed.
3. Probability-current decomposition
   `j/rho = F - D grad(log rho)`, an EPR-like circulation magnitude, current
   divergence residual, and an independent coarse Markov detailed-balance test.
   The continuous-current branch is QC-only and cannot determine the formalism
   when its stationarity-divergence check fails.
4. Two different Fisher quantities:
   - state-space Fisher/score energy;
   - path Fisher metric for `xyz`, `xy`, `z`, and the conditional slow-variable
     contribution.
5. Centered Fisher geometry and dynamic friction. The correlation time comes
   from the trajectory, while covariance amplitude comes from the same density
   and centered score that define path FI.
6. Exact finite-state and sampled Hatano–Sasa checks plus a sampled forward/reverse
   path-probability IFT. These are labelled as NESS relations, not as physical
   work identities.
7. A strict formalism verdict. Classical Jarzynski–Crooks remains blocked unless
   detailed balance and a physical beta/energy/conjugate-work mapping all pass.

The full frozen design is 264 supported scenarios, 32 dependent cell-pair
combinations, 31 path points, and five seeds: 40,920 stationary SDE runs.
Version 1.0.1 additionally resolves these cells to four WT and two SCA3 animals
in the core-secure cohort and writes an eight-animal-pair sensitivity summary.

## Server installation

```bash
unzip NeuroThermo_nonequilibrium_geometry_v1_0_1.zip
cd NeuroThermo_nonequilibrium_geometry_v1_0_1
python3 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the smoke test first:

```bash
./smoke_test.sh
```

The publication release contains the analysis-ready inputs at
`../../../data/inputs`, so no machine-specific input path is required.

```bash
./start_nohup.sh
tail -f nonequilibrium_geometry_v1_0_1.log
```

When the run is complete:

```bash
./pack_results.sh
```

Upload `neurothermo_nonequilibrium_geometry_results_v1_0_1.zip`. The archive
excludes resumable JSON checkpoints but includes compressed Markov caches,
tables, figures, manifests, and the formalism verdict.

## Interpretation rule

`FORMALISM_VERDICT.json` is authoritative. If it says `NESS`, then `phi` is a
nonequilibrium potential and **not** a physical Hamiltonian. A successful exact
Hatano–Sasa identity is a code/Markov consistency check; it does not by itself
validate sampled exponential averages or Jarzynski–Crooks.

## Corrections relative to v1.0.0

- All formalism decisions use the Markov time-reversal test; the failed
  continuous-current estimator is labelled diagnostic-invalid.
- Adaptive and linear protocols use 15 strictly unique positions each.
- Cycle affinities retain scenario, pair, and seed provenance.
- Cell-pair and animal-pair balanced outputs are both written.
- Endpoint membership is explicit, allowing exclusion of non-oscillatory
  endpoint families without guessing from representative IDs.
