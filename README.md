# AnnotaPipeline

AnnotaPipeline annotates hypothetical eukaryotic proteins by combining gene
prediction, sequence similarity and functional evidence. Transcriptomics and
proteomics can be added when RNA-seq or MS/MS data are available.

The project was developed by the Laboratório de Bioinformática at the
Universidade Federal de Santa Catarina, Brazil. It was tested on Unix-like
systems; Windows is not supported or tested.

There are two entry points:

- a Nextflow DSL2 workflow for scheduling, resume and published outputs;
- the legacy `scripts/annotapipeline.py` runner for direct execution.

Both use the same Python parsers and classification helpers, so the scientific
rules are identical.

## Contents

- [What changed](#what-changed)
- [Pipeline stages](#pipeline-stages)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Running the pipeline](#running-the-pipeline)
- [Outputs](#outputs)
- [Validation status](#validation-status)
- [Historical commits](#historical-commits)
- [Repository layout](#repository-layout)
- [Development](#development)
- [Citation and contact](#citation-and-contact)

## What changed

The v1 baseline is commit `b8bfc65`. The main changes are:

- **Configuration:** validated TOML loaded with Python 3.12 `tomllib` replaced
  YAML application files with path placeholders.
- **Python style:** Typer CLIs, `pathlib`, `Bio.SeqIO`, f-strings and shared
  helpers replaced mixed string handling and repeated setup.
- **Orchestration:** a Nextflow DSL2 workflow was added; the direct runner was
  retained.
- **Logging:** one shared module replaced repeated per-script setup.
- **Search:** BLAST+ remains the default; DIAMOND is optional for SwissProt
  and secondary searches.
- **Tabular output:** `.outfmt6` names became `.tsv`, with the same
  nine-column contract.
- **Environment:** Pixi, Ruff, `ty`, CI and an orchestrator-only Docker image
  replaced manual or Conda-only setup.
- **Tests:** 35 fixture tests now cover the CLI, configuration, parsers,
  runner and workflow.

Scientific behaviour was preserved:

- protein, genome and protein-plus-GFF input modes;
- SwissProt search followed by one secondary database;
- BLAST coverage `>`, identity and positivity `>=`; ties keep the first
  highest-bitscore hit;
- Kallisto TPM `>=` and Percolator q-value `<=`;
- GFF3/AUGUSTUS handling and the published output schemas.

Native YAML files remain only where their tools require them:
`annota_environment.yaml` for Conda and `.github/workflows/workflow-tests.yml`
for GitHub Actions.

DIAMOND is optional. It needs a prebuilt `.dmnd` database; no automatic
`diamond makedb` step is performed. DIAMOND has no `qcovs` field, so coverage
uses `qcovhsp`. BLAST+ runs are unchanged. RPS-BLAST/CDD remains BLAST+-only.

## Pipeline stages

- **Gene prediction:** AUGUSTUS predicts genes from a genome. Proteins, coding
  sequences and exons are extracted from GFF3 with `getAnnoFasta.pl`.
- **Similarity:** cleaned proteins are searched against SwissProt, then one
  secondary database. Hits are filtered by keywords, identity, positivity
  and coverage.
- **Functional annotation:** InterProScan runs on annotated and hypothetical
  proteins. HMMER/Pfam and RPS-BLAST/CDD add evidence for hypothetical
  proteins.
- **Transcriptomics:** optional Kallisto quantification of predicted coding
  sequences; TPM filtering uses the mean, median or a fixed value.
- **Proteomics:** optional Comet MS/MS and Percolator processing; PSMs are
  filtered by q-value and aggregated per protein.
- **Summary:** annotation sources are joined into product tables, annotated
  FASTA/GFF3 and a per-protein summary table.

Detailed stage descriptions are in [`docs/scripts.md`](docs/scripts.md).

## Requirements

Pixi provides the pinned runtime: Python 3.12.16, Nextflow 25.10.4, OpenJDK
21, Biopython 1.85, pandas 2.2.3 and Typer 0.29.2.

External scientific executables are still required:

- BLAST+ and RPS-BLAST;
- DIAMOND, only when `aligner = "diamond"`;
- AUGUSTUS and `getAnnoFasta.pl`, for genome input;
- InterProScan and HMMER;
- Kallisto, for transcriptomics;
- Comet and Percolator, for proteomics.

Required reference data are SwissProt, Pfam and CDD. Choose TrEMBL, EuPathDB,
NR or a custom FASTA database as the secondary database. DIAMOND databases
must already be formatted. Nextflow reference directories must contain a
non-empty `SHA256SUMS` manifest; see the
[reference contract](docs/nextflow.md#reference-contract).

The Docker image supplies only Python, Nextflow, Java and Python libraries.
Scientific tools and databases are excluded and must be mounted under `/data`
or configured per site.

## Quick start

```bash
git clone https://github.com/labioinfoufsc/AnnotaPipeline.git
cd AnnotaPipeline
pixi install --locked
cp config/config_example.toml my-run.toml
```

Edit `my-run.toml` before running. `config/annotapipeline.toml` is a template
with placeholders and intentionally fails validation until paths, organism
and basename are replaced. `config/config_example.toml` is a worked example.

Check the environment:

```bash
pixi run --locked python scripts/annotapipeline.py -h
```

TOML is the only application configuration format. Old YAML application files
are rejected; [`config/migration.md`](config/migration.md) describes a one-off
conversion. With `threads = 0`, the runner reads `SLURM_CPUS_PER_TASK` or
`NSLOTS`, then falls back to one thread.

## Running the pipeline

### Nextflow workflow

```bash
pixi run --locked nextflow run . \
  -params-file params.json \
  -output-dir results \
  -resume
```

`params.json` is a native Nextflow JSON parameter file; Nextflow does not read
TOML. Exactly one of `genome` or `protein` is required, and `gff` is valid
only with `protein`. Published outputs are written to `results/<basename>/`.

Resource defaults are conservative. Set `cpus`, `memory` and `time` in a site
config or through labels. The `docker`, `apptainer` and `conda` profiles
enable execution engines; they do not bundle scientific tools or databases.
See [`docs/nextflow.md`](docs/nextflow.md) for the full parameter reference.

### Legacy runner

Protein input:

```bash
pixi run --locked python scripts/annotapipeline.py \
  -c my-run.toml -p proteins.fasta
```

Genome input:

```bash
pixi run --locked python scripts/annotapipeline.py \
  -c my-run.toml -s genome.fasta
```

Protein plus GFF3 input:

```bash
pixi run --locked python scripts/annotapipeline.py \
  -c my-run.toml -p proteins.fasta -gff predictions.gff
```

The genome mode runs AUGUSTUS. The protein-plus-GFF mode skips prediction and
uses the supplied GFF3. Protein-only mode does not produce an annotated GFF.

## Outputs

Five main files are written to the launch directory:

- `All_Annotated_Products.txt`: identifiers and final annotations;
- `AnnotaPipeline_<basename>_proteins.fasta`: annotated proteins;
- `<basename>_Annotated_GFF.gff`: annotated GFF3, except in protein-only mode;
- `AnnotaPipeline_<basename>_transcripts.fasta`: annotated coding sequences;
- `AnnotaPipeline_<basename>_Summary.tsv`: per-protein summary table.

Raw outputs are kept in numbered folders:

- `1_GenePrediction_<basename>`;
- `2_SimilarityAnalysis_<basename>`;
- `3_FunctionalAnnotation_<basename>`;
- `4_TranscriptQuantification_<basename>`;
- `5_PeptideIdentification_<basename>`.

Output files for organisms in the publication are available at
<http://150.162.6.129/Annotafiles/>.

## Validation status

Validation was run with Pixi 0.63.2-compatible locked dependencies on Linux.

- **Format:** pass; `pixi run --locked format`, 22 files unchanged.
- **Lint and types:** pass; `pixi run --locked lint`, Ruff and `ty` 0.0.108.
- **Tests:** pass; `pixi run --locked test`, 35 fixture tests.
- **Nextflow config:** pass; `nextflow config -profile local` parses.
- **Container:** pass; rebuild plus Python, Nextflow, CLI help and config
  smoke checks.
- **Patch hygiene:** pass; `git diff --check` after the final edits.

The tests replace external tools with deterministic doubles. They cover CLI
contracts, configuration validation, parser invariants, all three direct-runner
input modes, workflow DAG behaviour, output publishing and resume.

Not yet validated:

- real BLAST+/DIAMOND/AUGUSTUS/HMMER/InterProScan/Kallisto/Comet/Percolator
  runs against reference databases;
- byte-level equivalence with v1 outputs on real biological data;
- Windows execution.

For local `ty` runs, a user-level `~/.config/ty/ty.toml` that injects a missing
`src` path must be isolated, for example with `XDG_CONFIG_HOME`. CI and Docker
are unaffected.

## Historical commits

Use `git show <hash>` to inspect each modernisation step.

- `b8bfc65` (2023-02-09): v1 baseline (`upstream/v1.0`), legacy runner and
  YAML configuration.
- `2e4b7c0` (2026-11-23): v2.0.0, Nextflow DSL2 workflow, lowercase scripts,
  shared helpers, Pixi, docs and tests.
- `89c5dad` (2026-11-23): modernisation, TOML, Typer, DIAMOND, shared logging,
  `.tsv` outputs, `ty`, Docker and expanded tests.

## Repository layout

```text
scripts/       Python CLIs, Nextflow task adapters and shared helpers
config/        TOML templates and migration guide
lib/           Nextflow parameter and metadata helpers
modules/       Nextflow processes
subworkflows/  Optional omics subworkflow
workflows/     Main DSL2 workflow
tests/         Fixture-based test suite
docs/          Detailed script and Nextflow references
```

## Development

```bash
pixi run --locked format
pixi run --locked lint
pixi run --locked test
```

Ruff uses the settings in `ruff.toml`, including an 88-character line length.
`ty` uses `ty.toml` so the flat `scripts/` module layout resolves correctly.
Comments use lowercase, explanatory prose and avoid decorative separators or
text that simply restates the code.

Command flags, inputs and outputs are listed in
[`docs/scripts.md`](docs/scripts.md).

## Citation and contact

If AnnotaPipeline was useful in your research, cite:

> AnnotaPipeline: An integrated tool to annotate eukaryotic proteins using
> multi-omics data. <https://doi.org/10.3389/fgene.2022.1020100>

Please also cite the external tools used in your run.

Questions and suggestions: <labioinfo.genome@gmail.com> or through GitHub
issues.

Licensed under the Apache License 2.0.
