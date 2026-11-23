# Configuration

The legacy Python runner uses `annotapipeline.toml` as its only configuration
format. Copy it, replace the path and organism placeholders, then pass it with
`-c`:

```bash
pixi run --locked python scripts/annotapipeline.py \
  -c config/annotapipeline.toml -p proteins.fasta
```

`config_example.toml` is a worked example. `annota_environment.yaml` is a Conda
convenience environment; the Pixi workspace in `pixi.toml` is the pinned
runtime.

The `[local_aligner]` section selects `aligner = "blast"` or `"diamond"`.
`blast-exe` is the executable for the selected aligner. DIAMOND searches need
DIAMOND-formatted reference databases and use `qcovhsp` for the coverage
threshold because DIAMOND does not emit BLAST's `qcovs` field.

`threads = 0` reads `SLURM_CPUS_PER_TASK` or `NSLOTS`, then falls back to one
thread. A positive integer overrides scheduler detection.

## Nextflow

The Nextflow workflow does not use these files. It uses a native JSON
`-params-file`; see [`docs/nextflow.md`](../docs/nextflow.md).
