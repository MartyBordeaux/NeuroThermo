# Validation targets

The package is expected to reproduce:

- 18 primary multi-sweep cells: 12 WT and 6 SCA3;
- q=.75 experimental support for all 18 cells;
- q=.50 support for 17 cells, with SCA3_06 absent rather than imputed;
- equal total biological weight per cell within genotype after splitting weight across near-optimal solutions;
- group support weights summing to 1 separately for WT and SCA3;
- no new parameter fitting and no formal animal-level genotype p-values.

## Packaging/import contract (v1.0.1)

`endpoint_v1/` is located at the repository root. Therefore the documented commands

```bash
python3 -m endpoint_v1.cli validate --config configs/server_endpoint_v1_0.yaml
./run_server.sh configs/server_endpoint_v1_0.yaml
```

work directly from the unpacked top-level directory without `pip install -e .`.
`run_server.sh` also changes to its own directory before execution.
