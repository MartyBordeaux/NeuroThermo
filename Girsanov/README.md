# Girsanov / NeuroThermo (2026-10-12)

This directory archives the finite-time Girsanov analysis of the **available** WT–SCA3 Hindmarsh–Rose retained-parameter ensemble. It is an extension, **not** a modification of the frozen publication benchmark.

- 8 WT and 4 SCA3 core-secure cells; 32 dependent endpoint pairs; 264 retained parameter scenarios.
- 11 prescribed analytical routes on an 11x11 intrinsic/drive grid; 525 frozen-route analysis was **not** reproduced because its binary route array was not available.
- Direct finite-horizon endpoint path KL for matched-current and zero-current controls; noise conditions 0.5D, D, 2D.
- Model-only information differences; these are **not** physical heat, entropy production, or a calibrated biological marker.
- The fitted slow component z is generally not dominant and the path attribution is parameter-nonidentifiable.

Scientific report: [`SCIENTIFIC_REPORT_RU.md`](SCIENTIFIC_REPORT_RU.md). Detailed methods and inventory: [`README_RU.md`](README_RU.md).

Source archive from the calculation session: `NeuroThermo_Girsanov_all_available_2026-10-12.zip` (49 files, SHA-256 `6ad261e543a90aa08e399dc3e1279bef6c78e36ab86f5742466d0ef4ad7815fd`). **The archive is not yet stored in GitHub**; individual available outputs are being added here. Missing files must not be inferred from this manifest.

Never count pairs, retained fits, routes, or stochastic trajectories as independent animals. No refit, refolding, or retuning was performed.
