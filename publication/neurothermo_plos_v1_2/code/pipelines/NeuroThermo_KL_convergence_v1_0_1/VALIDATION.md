# Validation

The package validates the explicit frozen directory and refuses automatic
input discovery. The server configuration requires exactly 264 selected
scenarios, 32 dependent pairs, 31 path positions, the frozen five-seed
sequence, the three declared time steps, and the original stationary-density
settings.

The server configuration additionally requires the `all_task_extrema` grid
strategy, a positive coverage margin, and exact 100% retained mass. The smoke
test performs both deterministic passes and asserts
`all_tasks_full_coverage=true` and `minimum_retained_mass=1.0`.

Checkpoints contain the scientific fingerprint and are rejected after any
change to the scientific configuration, frozen input files, animal mapping, or
package version. Worker count, output directory, and resume setting do not
change the fingerprint.

Pilot and analysis checkpoints are separate but share the same scientific
fingerprint. The smoke calculation verifies the complete two-pass data path,
coverage audit, verdict tiers, and output schema using synthetic inputs. Its
scientific verdict is disabled and must not be cited.
