# Configuration migration: YAML to TOML

`config/annotapipeline.toml` is the only configuration accepted by
`annotapipeline.py`. Python 3.12 parses it with `tomllib`, so no YAML runtime
dependency is required. TOML preserves types explicitly: `on`/`off` stay
strings instead of becoming YAML booleans, and integers, floats and lists
cannot be confused with text. This removes a real defect in the legacy runner,
where an `on` value was rendered as `True` in AUGUSTUS command lines.

The Nextflow workflow continues to use a native JSON `-params-file`. Its
input validation, reference staging and cache invalidation are separate
contracts; converting that interface to TOML would add a custom parser with
no benefit. GitHub Actions, `annota_environment.yaml` and Nextflow's own
configuration files remain YAML/native Groovy.

## Convert a legacy YAML file once

The loader no longer accepts YAML. Install conversion tools temporarily,
write a new `.toml` file and keep the original YAML as an external archive:

```bash
python -m pip install 'tomli-w>=1,<2' PyYAML
python - <<'PY'
import tomli_w
import yaml
from pathlib import Path

source = Path("archived/annotapipeline.yaml")
target = Path("config/annotapipeline.toml")
data = yaml.safe_load(source.read_text())
renamed = {
    "augustus-optional": "augustus_optional",
    "seq-cleaner": "seq_cleaner",
    "local-aligner": "local_aligner",
}
data = {renamed.get(key, key): value for key, value in data.items()}
# Replace nulls with TOML-compatible empty strings, lists or zeros as
# appropriate before writing.
target.write_text(tomli_w.dumps(data))
PY
```

A mechanical rename is a starting point, not sufficient on its own. Check
each of these semantic differences:

| YAML input | TOML replacement |
| --- | --- |
| `value:` or `value: null` | `value = ""` for text, `[]` for lists, or omit the key |
| `on` / `off` unquoted | `on = "on"` / `on = "off"` (quotes are mandatory) |
| `E: 1e-5` | `E = 1e-5` (TOML floating point) |
| `rna-seq: [file1, file2]` | `rna-seq = ["file1", "file2"]` |
| `augustus-optional:` | `[augustus_optional]` |
| `seq-cleaner:` | `[seq_cleaner]` |
| `local-aligner:` | `[local_aligner]` |
| `first:` / `last:` left blank | omit both keys from `[proteomics]` |
| `ms:`, `l:` or `s:` left blank | omit the key from its section |

Do not translate blank numeric settings to zero: `first = 0` and
`last = 0` would override the scan range in `comet.params`, while `l = 0`
and `s = 0` would be passed to Kallisto.

`threads = 0` reads `SLURM_CPUS_PER_TASK` or `NSLOTS`, falling back to one
thread. A positive integer overrides this behaviour.

## New settings

Add the aligner to `[local_aligner]`:

```toml
[local_aligner]
aligner = "blast"  # or "diamond"
blast-exe = "blastp"
```

DIAMOND databases must be prepared with `diamond makedb`; `blast-exe` then
points to the `diamond` executable. The coverage threshold uses DIAMOND's
`qcovhsp` field because `qcovs` is not available.

Tabular search outputs now use `.tsv` names instead of `.outfmt6`. The column
order and scientific thresholds are unchanged for BLAST+ runs. Output
directory names, file schemas and scientific steps depend on setting values,
not on the file format.
