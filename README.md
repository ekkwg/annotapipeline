[![Quote](<https://img.shields.io/badge/Quote-May%20the%20Dark%20shine%20your%20way%20(Darkdiver%20Grandahl)-blueviolet>)](https://g.co/kgs/bMJLfj)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

<font size=16>**AnnotaPipeline**</font>

- [**About AnnotaPipeline**](#about-annotapipeline)
- [**Analysis description**](#analysis-description)
  - [**1. Gene prediction**](#1-gene-prediction)
  - [**2. Similarity analysis**](#2-similarity-analysis)
  - [**3. Functional annotation**](#3-functional-annotation)
  - [**4. Optional: Transcript quantification**](#4-optional-transcript-quantification)
  - [**5. Optional: Peptide identification**](#5-optional-peptide-identification)
  - [**Last but not least**](#last-but-not-least)
- [**How to setup**](#how-to-setup)
  - [**Download databases**](#download-databases)
  - [**InterProScan: Tested databases**](#interproscan-tested-databases)
  - [**Using Pixi**](#using-pixi)
  - [**Python scripts**](#python-scripts)
  - [**Using conda**](#using-conda)
    - [**Setup AUGUSTUS species for personalized predictions**](#setup-augustus-species-for-personalized-predictions)
  - [**Manual install**](#manual-install)
  - [**Configure the runner**](#configure-the-runner)
  - [**Container image**](#container-image)
- [**How to run**](#how-to-run)
  - [**Nextflow workflow**](#nextflow-workflow)
  - [**Python runner**](#python-runner)
    - [**Protein file as input**](#protein-file-as-input)
    - [**Genome file as input**](#genome-file-as-input)
    - [**Protein and GFF files as input**](#protein-and-gff-files-as-input)
- [**Output**](#output)
- [**Output Files**](#output-files)
- [**Tests**](#tests)
- [**Citation**](#citation)


# **About AnnotaPipeline**

Integrated tool to annotate hypothetical proteins developed by Laboratório de Bioinformática at Universidade Federal de Santa Catarina (Brazil).

AnnotaPipeline was tested in Unix-based systems. We strongly recommend to use Unix-based systems or Mac. AnnotaPipeline was not tested in Windows.

If you have questions, suggestions or difficulties regarding the pipeline, please do not hesitate to contact our team here on GitHub or by email: <labioinfo.genome@gmail.com>.

# **Analysis description**

> `annotapipeline.py` checks all required parameters in `config/annotapipeline.toml` before execution

## **1. Gene prediction**

Run AUGUSTUS for gene prediction, using these required arguments:

- strand
- genemodel
- species
- protein
- introns
- start
- stop
- cds

You can also use these optional arguments:

> Or add AUGUSTUS arguments below `[augustus_optional]` (as described in the example config)

- hintsfile
- extrinsicCfgFile
- UTR

After gene prediction, sequences are "cleaned" using `minsize-seq` from `[seq_cleaner]` in `config/annotapipeline.toml`.

`.aa`, `.cdsexon` and `.codingseq` sequences are extracted from GFF output using `getAnnoFasta.pl` (AUGUSTUS script).

`.aa` sequences are used for subsequent analysis. `.codingseq` sequences are used for transcriptomics analysis (optional).

## **2. Similarity analysis**

Similarity analysis runs through `blastp_parser.py`, which executes BLAST+ or DIAMOND with tab-separated output (parsed after the run) on cleaned predicted proteins against SwissProt:

- qseqid
- sseqid
- sacc
- bitscore
- evalue
- ppos
- pident
- qcovs
- stitle
  > `evalue` and `max_target_seqs` given by user

This output is parsed to find annotations in the secondary database. The keyword list in `[pipeline]` is used to exclude potential hypothetical annotations.

Hits are classified &ndash; as a potential annotation &ndash; if: (i) it doesn't have any words in the keyword list, and (ii) passed value thresholds for identity, positivity and coverage. All potential annotations &ndash; hits that passed all criteria &ndash; are present in `BASENAME_SwissProt_annotations.txt` for manual check.

Potential annotations are ranked by `bitscore` and the best hit is assigned to the corresponding protein. Proteins that aren't annotated based on SwissProt annotations, are separated in `BASENAME_BLASTp_AA_SwissProted.fasta` and the aligner is rerun against the secondary database (selected in `[databases]`).

Set `[local_aligner] aligner = "blast"` (default) or `"diamond"`. For DIAMOND, `blast-exe` points to the `diamond` binary and `swissprot-db`/`secondary-db` must be DIAMOND-formatted databases. DIAMOND has no `qcovs` field, so the coverage threshold uses its `qcovhsp` output; BLAST+ runs are unchanged.

Following the same process, `BASENAME_SpecifiedDB_annotations.txt` classifies potential annotations.

Annotation files:

- `BASENAME_annotated_products.txt` contains all annotations
- `BASENAME_hypothetical_products.txt` contains hypothetical proteins with no hits that passed all criteria(which doesn't have a single hit that passes all criteria) are present in `BASENAME_hypothetical_products.txt`
- `BASENAME_no_hit_products.txt` contains proteins with no hits &ndash; against SwissProt and the secondary database &ndash; that are treated as hypothetical for subsequent analysis

## **3. Functional annotation**

**Functional analysis runs in two different ways:**

- For annotated proteins:
  - InterProScan is used (with configured databases) to get ontology and IPR terms
- For hypothetical proteins (which includes no_hit_products):
  - InterProScan, RPS-BLAST and HMMER are used to find hits of possible functions of predicted proteins

**Software arguments:**

> Optional arguments can be given in the config file for InterProScan, HMMER and RPS-BLAST (not tested)

- InterProScan &ndash; for annotated and hypothetical proteins &ndash; uses `-goterms` and `-iprlookup` arguments.
- `hmmscan` runs with `--noali` argument and user values for `evalue` and `domE`. It also uses Pfam database.
- RPS-BLAST runs with `evalue` and `max_target_seqs` from `[rpsblast]` and the same nine-column tabular contract written to `BASENAME_rpsblast_output.tsv`. RPS-BLAST uses the CDD database and remains BLAST+-only.

**Parsing:**

- `funcannotation_parser.py` join outputs for both InterProScan runs in a single file called `InterProScan_Out_BASENAME.txt`. This file summarizes outputs for each predicted protein.
  - Coils, Gene3D and MobiDBLite databases are structural databases and are excluded from this output.
- RPS-BLAST information for hypothetical proteins are summarized in `BASENAME_Grouped_Hypothetical_Information.txt` as it gives long descriptions. It may be helpful to find functional hints for proteins.
- `info_parser.py` uses InterProScan results with parsed files from BLAST results to generate `All_annotated_products.txt`. This file joins gene annotation (from BLAST) with functional annotation (from InterProScan) using GOs and IPR values. This file is used to annotate sequences (aminoacid and nucleotide) in FASTA and GFF files.

## **4. Optional: Transcript quantification**

Transcript quantification uses Kallisto and transcripts (`.codingseq` provided by AUGUSTUS) with RNA-seq data given by user.

Analysis starts with `kallisto index` and is followed by `kallisto quant`:

> Both methods uses `bootstrap` value provided by user

- 1 &ndash; paired-end data &ndash; using estimated average fragment length (`l`) and estimated standard deviation of fragment length (`s`)

  - Optional arguments for paired-end data in `[kallisto]`, if given

  OR

- 2 &ndash; single-end data &ndash; using **REQUIRED** `l` and `s` values from `[kallisto]`

**Parsing:**

- `kallisto_parser.py` removes hits &ndash; by TPM value &ndash; below the `threshold` parameter from `[kallisto]`. Possible thresholds are:
  - Median
  - Mean
  - Float value (user input)
- Parsed output file is called `BASENAME_Transcript_Quantification.tsv` and contains a simplified result of Kallisto output (`abundance.tsv`) with target_id and TPM value

## **5. Optional: Peptide identification**

> **WARNING 1:** Failed Comet MS/MS runs do not stop AnnotaPipeline, so check `AnnotaPipeline_LogFile.log` to confirm that every spectrometry file produced an output

> **WARNING 2:** Before run, check if your comet.params file is compatible with installed comet version

Proteomics analysis uses Comet MS/MS with `comet.params` config given by user. In this file, our script overwrite values for following parameters:

- `decoy_search = 1`
- `output_pepxmlfile = 0`
- `output_percolatorfile = 1`
- `decoy_prefix = DECOY_`

Comet MS/MS runs with:

- _Modified_ `comet.params`
- Annotated protein file
- Path containing mass spectrometry files (`comet-spectrometry` param) and extension (`comet-ext` param)
  - Extension can be _mzXML, mzML, Thermo raw, mgf, and ms2 variants (cms2, bms2, ms2)_
- Optional arguments `first` and `last` are used, if given, and overwrite cutoff values defined in `comet.params`

Comet MS/MS outputs (files with `.pin` extension) and raw files are moved to current directory.

Percolator runs for each `.pin` file with default parameters.

Percolator outputs are parsed by `percolator_parser.py` using `percolator-qvalue` from `[proteomics]`. Raw files are kept. `quantitative_proteomics` uses parsed Percolator files to create `BASENAME_Total_Proteomics_Quantification.tsv`.

This output quantifies (all spectometry outputs):

- Unique Peptide &ndash; number of unique peptides found across the entire dataset
- Total Peptide &ndash; number of total peptides found across the entire dataset
- Unique Spectrum &ndash; number of unique spectrum found across the entire dataset
- Total Spectrum &ndash; number of unique spectrum found across the entire dataset

## **Last but not least**

`summary_parser.py` integrates all outputs (in a single file summarizing all annotations found for each protein) from:

- Prediction and similarity analysis
- Functional annotation
- Transcriptomic (if used)
- Peptide identification (if used)

# **How to setup**

AnnotaPipeline requires the following software to run properly:

- BLAST+ and RPS-BLAST (available at <https://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST>)
- DIAMOND (optional, <https://github.com/bbuchfink/diamond>) when `aligner = "diamond"`
- InterProScan (available at <https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html>)
- HMMER (available at <http://hmmer.org/download.html>)

You will also need to install [AUGUSTUS](https://doi.org/10.1093/bioinformatics/btr010) (available at <https://github.com/Gaius-Augustus/Augustus>) if you want to run this pipeline starting with gene/protein prediction.

Before executing, copy `config/annotapipeline.toml` and replace the path and
organism placeholders.

Repository directories are lower case: pipeline scripts live in `scripts/`
and configuration files in `config/`.

## **Download databases**

**Required:**

- SwissProt (available at <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz>)
- Pfam (available at <http://pfam.xfam.org>)
- CDD (available at <https://www.ncbi.nlm.nih.gov/cdd>)

**Choose one secondary database:**

- TrEMBL (available at <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz>)

- EuPathDB (available at <https://veupathdb.org/veupathdb/app>):

  - AmoebaDB
  - CryptoDB
  - FungiDB
  - GiardiaDB
  - HostDB
  - MicrosporidiaDB
  - PiroplasmaDB
  - PlasmoDB
  - ToxoDB (tested)
  - TrichDB
  - TriTrypDB (tested)

- NCBI | NR Database (available at <https://ftp.ncbi.nlm.nih.gov/blast/db>)

  > **TIP:** You can use a subset of the NR Database

- Custom Database: You can set a custom database if you provide the pattern to get descriptions in the fasta file.

  Example:

  Sequences from the Arabidopsis database ([ArabdopsisDB Tair10](https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_pep_20101214)) are separated by pipes ("`|`") with protein descriptions in the 3th field.

  ```
  >AT1G51370.2 | Symbols:  | F-box/RNI-like/FBD-like domains-containing protein | chr1:19045615-19046748 FORWARD LENGTH=346
  ```

  It is possible to use AnnotaPipeline with the Tair10DB changing the following parameters:

  - `secondary-format` set to `custom`
  - `customsep` set to `"|"`
  - `customcolumn` set to `2`

  > **NOTE:** `customcolumn` is set to `2` because indexing starts at 0.

## **InterProScan: Tested databases**

> WARNING: Installation through conda/mamba requires manual download and configuration of InterProScan databases

- CDD
- Gene3D
- Hamap
- Panther
- Pfam
- Pirsf
- Pirsr
- Prints
- PrositePatterns
- PrositeProfiles
- Sfld
- Smart (unlicenced)
- Superfamily
- Tigrfam
- Coils
- MobiDBLite

## **Using Pixi**

[Pixi](https://pixi.sh) is the pinned runtime for both entry points. It
provides Python `3.12.16`, Nextflow `25.10.4`, OpenJDK `21.0.13.1`, Biopython
`1.85`, pandas `2.2.3` and Typer `0.29.2`. The development environment adds
Ruff and `ty`. Install Pixi, then run:

```bash
pixi install --locked
pixi run --locked python scripts/annotapipeline.py -h
```

The `runtime` Pixi environment omits Ruff and `ty` and is used by the Docker image.
External tools and reference databases are still required. Each Nextflow
reference directory must contain a non-empty `SHA256SUMS` manifest; see the
[reference contract](docs/nextflow.md#reference-contract).

## **Python scripts**

The Python layer in `scripts/` contains the Nextflow task adapters, shared
I/O and configuration helpers, and the command-line parsers. Each CLI accepts
`-h`, `-help`, and `--help`. Identifiers follow Python `snake_case` and file
names are lower case. Ruff handles formatting and style; `ty` checks types:

```bash
pixi run --locked format
pixi run --locked lint
```

Command flags, inputs, and outputs are listed in
[`docs/scripts.md`](docs/scripts.md).

## **Using conda**

The Conda environment is a convenience option for the direct runner; Pixi is
the pinned runtime. From the repository root:

```bash
conda env create -n annotapipeline -f annota_environment.yaml
conda activate annotapipeline
```

[mamba](https://github.com/mamba-org/mamba) can replace `conda env create`
with `mamba env create`. Then [**download databases**](#download-databases),
configure [InterProScan](#interproscan-tested-databases), and follow
[configure the runner](#configure-the-runner).

### **Setup AUGUSTUS species for personalized predictions**

1. Locate the environment home:

```bash
echo "$CONDA_PREFIX"
```

2. Add the trained species folder under `$CONDA_PREFIX/config/species`.

## **Manual install**

1. Clone the repository:

```bash
git clone https://github.com/labioinfoufsc/AnnotaPipeline.git
cd AnnotaPipeline
```

2. Install the Python package from the repository root; scripts in `scripts/`
   become available on `$PATH`:

```bash
python3 -m pip install .
```

3. Install the required software:

   - BLAST+ and RPS-BLAST (available at <https://ftp.ncbi.nih.gov/blast/executables/blast+/LATEST>)
- DIAMOND (optional, <https://github.com/bbuchfink/diamond>) when `aligner = "diamond"`
   - InterProScan (available at <https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html>)
   - HMMER (available at <http://hmmer.org/download.html>)

4. Optional software:

   - Kallisto (available at <https://pachterlab.github.io/kallisto/download.html>)
   - Comet MS/MS (available at <https://github.com/UWPR/Comet/releases/latest>)
     - Requires Percolator (available at <https://github.com/percolator/percolator>)

5. [**Download databases**](#download-databases)

6. Follow [configure the runner](#configure-the-runner).

## **Configure the runner**

`config/annotapipeline.toml` is the canonical configuration for
`annotapipeline.py`. `config/config_example.toml` is a worked example. TOML is
parsed with Python 3.12's `tomllib`, so `on` and `off` stay strings and
numeric settings keep their types.

YAML application configuration is no longer supported; use a `.toml` file.
See [migration.md](config/migration.md) for a one-off YAML-to-TOML conversion.
The conversion is non-destructive: keep the original YAML file outside the
runner and write a new `.toml` file. GitHub Actions and
`annota_environment.yaml` remain YAML because their tools require it.

`threads = 0` reads `SLURM_CPUS_PER_TASK` or `NSLOTS`, falling back to one
thread. A positive integer overrides scheduler detection.

> **TIP:** InterProScan can also be installed separately and configured with
> `interproscan-exe`.

## **Container image**

`Dockerfile` builds an orchestrator-only image with the locked Pixi `runtime`
environment: Python, Nextflow, Java and Python libraries. BLAST+, DIAMOND,
AUGUSTUS, HMMER, InterProScan, Kallisto, Comet, Percolator and reference
databases are **not** included; mount them under `/data` or supply them through
site configuration.

```bash
docker build -t annotapipeline:orchestrator .
docker run --rm annotapipeline:orchestrator python --version
docker run --rm annotapipeline:orchestrator nextflow -version
```

The container runs as UID/GID 10001 with `/work` as its working directory.
Bind-mounted paths must be writable by that identity, or run with
`--user "$(id -u):$(id -g)"`. On SELinux hosts append `:z` or `:Z` to bind
mounts; otherwise the daemon can deny reads even to container root:

```bash
docker run --rm -v "$PWD/tests:/opt/annotapipeline/tests:ro,z" \
  -w /opt/annotapipeline annotapipeline:orchestrator \
  python -m unittest discover -s tests
```

# **How to run**

## **Nextflow workflow**

The Nextflow DSL2 workflow schedules tasks and manages `-resume` and
concurrency:

```bash
pixi run --locked nextflow run . \
  -params-file params.json \
  -output-dir results \
  -resume
```

`params.json` is a native Nextflow JSON parameter file; paths resolve relative
to the launch directory. Exactly one of `genome` or `protein` is required, and
`gff` is valid only with `protein`. Published outputs are written to
`results/<basename>/`. Resource defaults (`cpus`, `memory`, `time`) are
conservative; set them from the scheduler or a site config. The `docker`,
`apptainer` and `conda` profiles enable engines only, without binding process
containers or environments.

See [`docs/nextflow.md`](docs/nextflow.md) for the full parameter reference,
reference contract, execution profiles, and resume behaviour.

> **Validation:** fixture tools cover workflow behaviour and CLI contracts.
> Runs with the real toolchain, reference databases, and the v1.0 output
> baseline are pending.

## **Python runner**

The `annotapipeline.py` CLI reads `config/annotapipeline.toml` and accepts
three input modes. With Pixi:

- [**Protein file as input**](#protein-file-as-input)
- [**Genome file as input**](#genome-file-as-input)
- [**Protein and GFF files as input**](#protein-and-gff-files-as-input)

## **Protein file as input**

```bash
pixi run --locked python scripts/annotapipeline.py \
  -c config/annotapipeline.toml -p protein_sequences.fasta
```

The pipeline annotates the submitted `protein_sequences.fasta` and writes
FASTA headers with the organism and annotation description. This mode does
not produce an annotated GFF file.

## **Genome file as input**

```bash
pixi run --locked python scripts/annotapipeline.py \
  -c config/annotapipeline.toml -s genomic_data.fasta
```

AUGUSTUS predicts genes and proteins from `genomic_data.fasta`; the
predicted proteins enter similarity analysis. Use a trained AUGUSTUS model
for the species.

## **Protein and GFF files as input**

```bash
pixi run --locked python scripts/annotapipeline.py \
  -c config/annotapipeline.toml \
  -p protein_sequences.fasta -gff gff_file.gff
```

Supply `.aa` and `.gff` files from AUGUSTUS predictions. The submitted GFF
must use GFF3 format. Gene prediction is skipped; similarity analysis starts
from the supplied proteins.

# **Output**

AnnotaPipeline, it will output five main files (along with many others in their respective folders):

- `All_Annotated_Products.txt` contains all unique sequence identifiers and their respective annotations (with functional annotations &ndash; when present).
- `AnnotaPipeline_BASENAME_proteins.fasta` contains all sequences and their annotations (with functional annotations &ndash; when present) in FASTA format.
- `BASENAME_Annotated_GFF.gff` contains all sequences and their annotations (with functional annotations &ndash; when present) in GFF3 format.
  > This file is absent in [**Protein file as input**](#protein-file-as-input) run mode
- `AnnotaPipeline_BASENAME_transcripts.fasta` contains nucleotide sequences for predicted proteins, with the same features present in the protein file.
- `AnnotaPipeline_BASENAME_Summary.tsv` summarizes hits for each protein in similarity, functional, transcriptomics (if used) and proteomics analysis (if used).

Raw outputs are listed inside output folders:

- `1_GenePrediction_BASENAME` &ndash; AUGUSTUS files
- `2_SimilarityAnalysis_BASENAME` &ndash; BLAST+ or DIAMOND analysis
- `3_FunctionalAnnotation_BASENAME` &ndash; InterProScan/HMMER/RPS-BLAST analysis
- `4_TranscriptQuantification_BASENAME` &ndash; Kallisto analysis
- `5_PeptideIdentification_BASENAME` &ndash; Comet MS/MS and Percolator analysis

Output folders and files are written to the launch directory.

# **Output Files**

- Output files for Organisms present in AnnotaPipeline publication are available at: http://150.162.6.129/Annotafiles/


# **Tests**

```bash
pixi run --locked format
pixi run --locked lint
pixi run --locked test
```

Fixture tests replace external tools and cover configuration loading, parser
invariants, all three direct-runner input modes, DAG behaviour, output
publishing, resume, and CLI contracts. Real-tool validation against the real
databases and the v1.0 output baseline is still pending; see
[tests](docs/nextflow.md#tests).

# **Citation**

If you used AnnotaPipeline in your research, please cite us

[AnnotaPipeline: An integrated tool to annotate eukaryotic proteins using multi-omics data](https://doi.org/10.3389/fgene.2022.1020100)

...And the following papers:

- AUGUSTUS: [Stanke M. et al., 2003](https://doi.org/10.1093/bioinformatics/btg1080)
- BLAST+: [Camacho C. et al., 2008](https://doi.org/10.1186/1471-2105-10-421)
- HMMER: follow [HMMER user guide](http://eddylab.org/software/hmmer/Userguide.pdf)
- InterProScan: [Jones et al., 2014](https://doi.org/10.1093/bioinformatics/btu031)

If you used Transcriptomics module, please also cite:

- Kallisto: [Bray, N. L. et al., 2016](https://doi.org/10.1038/nbt.3519)

If you used Proteomics module, please also cite:

- COMET MS/MS: [Eng, J. K., et al., 2012](https://doi.org/10.1002/pmic.201200439)
- percolator: [The, M. et al., 2016](https://doi.org/10.1007/s13361-016-1460-7)
