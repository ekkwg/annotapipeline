#!/usr/bin/env python3
"""run the shell-driven AnnotaPipeline orchestrator"""

from __future__ import annotations

import contextlib
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import typer
from annotation_io import (
    annotate_codingseq,
    fasta_fetcher,
    quantitative_proteomics,
    sequence_cleaner,
)
from blast_classification import BLAST_FIELDS
from logging_utils import configure

from config import ConfigError, as_flag, load

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)
LOGGER = logging.getLogger("AnnotaPipeline")
ANNOTATIONS = "All_Annotated_Products.txt"
RAW_ANNOTATIONS = "All_annotation_products.txt"


class PipelineError(RuntimeError):
    """raise when a pipeline step cannot produce a required output"""


def configure_logging() -> None:
    configure("AnnotaPipeline_LogFile.log", name="AnnotaPipeline")


def run(command: list[str], stdout=None, *, allowed=(0,)) -> int:
    LOGGER.info(f"$ {' '.join(command)}")
    result = subprocess.run(command, stdout=stdout)
    if result.returncode not in allowed:
        raise PipelineError(
            f"command failed ({result.returncode}): {' '.join(command)}"
        )
    return result.returncode


def check_file(path: str | Path, *, empty: bool = False) -> None:
    target = Path(path)
    if not target.is_file():
        raise PipelineError(f"missing expected output: {target}")
    if not empty and target.stat().st_size == 0:
        LOGGER.warning(f"expected output is empty: {target}")


def script(name: str) -> Path:
    return Path(__file__).resolve().parent / name


def version_key(text: str) -> list[int | str]:
    """sort identifiers containing numbers in natural version order"""
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", text)
        if part
    ]


def sort_annotations() -> Path:
    lines = Path(RAW_ANNOTATIONS).read_text(encoding="utf-8").splitlines()
    target = Path(ANNOTATIONS)
    target.write_text(
        "".join(f"{line}\n" for line in sorted(lines, key=version_key)),
        encoding="utf-8",
    )
    Path(RAW_ANNOTATIONS).unlink(missing_ok=True)
    return target


def sort_tsv(source: Path, target: Path) -> Path:
    lines = source.read_text(encoding="utf-8").splitlines()
    if lines:
        header, rows = lines[0], lines[1:]
        content = "\n".join([header, *sorted(rows, key=version_key)]) + "\n"
    else:
        content = ""
    target.write_text(content, encoding="utf-8")
    source.unlink(missing_ok=True)
    return target


def flag_value(value: Any) -> str:
    """convert a legacy flag without changing free-form values"""
    text = str(value)
    return as_flag(value) if text.lower() in {"on", "off", "flag"} else text


def augustus_flags(settings: dict) -> list[str]:
    flags: list[str] = []
    for section in (settings["augustus"], settings["augustus_optional"]):
        for key, value in section.items():
            if key == "augustus-path" or value in (None, ""):
                continue
            if str(value).lower() == "flag":
                flags.append(f"-{key}")
            else:
                flags.append(f"--{key}={flag_value(value)}")
    return flags


def augustus_paths(settings: dict) -> tuple[Path, Path, Path]:
    configured = str(settings["augustus"]["augustus-path"])
    if configured.lower() == "conda":
        prefix = Path(os.environ.get("CONDA_PREFIX", ""))
        return (
            Path("augustus"),
            prefix / "config",
            Path("getAnnoFasta.pl"),
        )
    root = Path(configured)
    return (
        root / "bin" / "augustus",
        root / "config",
        root / "scripts" / "getAnnoFasta.pl",
    )


def tool_flags(settings: dict, *, skip: set[str]) -> list[str]:
    """translate optional tool settings into command-line flags"""
    flags: list[str] = []
    for key, value in settings.items():
        if key in skip or value in (None, ""):
            continue
        if str(value).lower() == "flag":
            flags.append(f"-{key}")
        else:
            flags.extend([f"-{key}", str(value)])
    return flags


def read_ids(path: str | Path) -> list[str]:
    return [
        line.split()[0]
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def interproscan(settings: dict, kind: str, basename: str) -> None:
    section = settings["interproscan"]
    output = f"{basename}_interproscan_{kind}_output.gff3"
    command = [
        str(section["interproscan-exe"]),
        "-i",
        f"{kind.capitalize()}_Products.fasta",
        "-o",
        output,
        "-f",
        "GFF3",
        "-t",
        "p",
        "-goterms",
        "-iprlookup",
        *tool_flags(section, skip={"interproscan-exe"}),
    ]
    LOGGER.info(f"Running InterProScan for {kind} proteins")
    run(command)
    check_file(output, empty=True)


def split_blast_products(
    basename: str, fasta: str | Path, similarity_dir: Path
) -> None:
    annotated = read_ids(similarity_dir / f"{basename}_annotated_products.txt")
    hypothetical = read_ids(similarity_dir / f"{basename}_hypothetical_products.txt")
    no_hit = read_ids(similarity_dir / f"{basename}_no_hit_products.txt")
    fasta_fetcher(fasta, hypothetical + no_hit, "Hypothetical_Products.fasta")
    fasta_fetcher(fasta, annotated, "Annotated_Products.fasta")
    for name in ("Annotated_Products.fasta", "Hypothetical_Products.fasta"):
        check_file(name)


def make_annotations(
    settings: dict,
    basename: str,
    similarity_dir: Path,
    functional_dir: Path,
    home: Path,
) -> None:
    run(
        [
            settings["pipeline"]["python"],
            str(script("funcannotation_parser.py")),
            "-ipr_annot",
            f"{basename}_interproscan_annotated_output.gff3",
            "-ipr_hyp",
            f"{basename}_interproscan_hypothetical_output.gff3",
            "-hmm",
            f"{basename}_hmmscan_output.txt",
            "-rpsblast",
            f"{basename}_rpsblast_output.tsv",
            "-basename",
            basename,
        ]
    )
    os.chdir(home)
    run(
        [
            settings["pipeline"]["python"],
            str(script("info_parser.py")),
            "-ipr1",
            str(functional_dir / f"{basename}_interproscan_annotated_output.gff3"),
            "-ipr2",
            str(functional_dir / f"{basename}_interproscan_hypothetical_output.gff3"),
            "-a",
            str(similarity_dir / f"{basename}_annotated_products.txt"),
            "-hy",
            str(similarity_dir / f"{basename}_hypothetical_products.txt"),
            "-nh",
            str(similarity_dir / f"{basename}_no_hit_products.txt"),
        ]
    )
    check_file(RAW_ANNOTATIONS, empty=True)
    sort_annotations()


def export_sequences(
    settings: dict,
    basename: str,
    protein: str | None,
    gff: str | None,
    augustus_gff: str | None,
    clean_fasta: str,
) -> None:
    python = settings["pipeline"]["python"]
    organism = settings["pipeline"]["organism"]
    if gff or (protein is None and augustus_gff):
        source_gff = gff or augustus_gff
        run(
            [
                python,
                str(script("gfftofasta_parser.py")),
                "-gff",
                str(source_gff),
                "-annot",
                ANNOTATIONS,
                "-b",
                basename,
                "-faf",
                clean_fasta,
                "-org",
                organism,
            ]
        )
        run(
            [
                python,
                str(script("fastatogff.py")),
                "-gff",
                str(source_gff),
                "-all",
                ANNOTATIONS,
                "-b",
                basename,
            ]
        )
        if protein is None and augustus_gff:
            annotate_codingseq(
                f"AnnotaPipeline_{basename}_proteins.fasta",
                Path(augustus_gff).with_suffix(".codingseq"),
                basename,
            )
    else:
        run(
            [
                python,
                str(script("fasta_simple.py")),
                "-annot",
                ANNOTATIONS,
                "-b",
                basename,
                "-faf",
                clean_fasta,
                "-org",
                organism,
            ]
        )
    check_file(f"AnnotaPipeline_{basename}_proteins.fasta")


def kallisto(
    settings: dict, basename: str, transcripts: Path | None, output_dir: Path
) -> bool:
    section = settings["kallisto"]
    if not section["kallisto-exe"]:
        return False
    if transcripts is None or not transcripts.is_file():
        LOGGER.warning("Kallisto skipped: no predicted coding sequences")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)
    executable = str(section["kallisto-exe"])
    reads = [str(Path(path).absolute()) for path in section["rna-seq"]]
    index = f"{basename}_kallisto_index.idx"
    out = Path(f"{basename}_kallisto_output")
    with contextlib.chdir(output_dir):
        run([executable, "index", "-i", index, str(transcripts)])
        check_file(index)
        command = [
            executable,
            "quant",
            "-i",
            index,
            "-o",
            str(out),
            "-b",
            str(section["bootstrap"]),
        ]
        if section["single_end"]:
            command.extend(
                [
                    "--single",
                    "-l",
                    str(section["l"]),
                    "-s",
                    str(section["s"]),
                ]
            )
        else:
            if section.get("l") not in (None, ""):
                command.extend(["-l", str(section["l"])])
            if section.get("s") not in (None, ""):
                command.extend(["-s", str(section["s"])])
        run([*command, *reads])
        abundance = out / "abundance.tsv"
        check_file(abundance)
        threshold = section["threshold"]
        if str(threshold).lower() == "mean":
            flags = ["-tpmavg"]
        elif str(threshold).lower() == "median":
            flags = ["-tpmmd"]
        else:
            flags = ["-tpmval", str(threshold)]
        run(
            [
                settings["pipeline"]["python"],
                str(script("kallisto_parser.py")),
                "-ktfile",
                str(abundance),
                "-basename",
                basename,
                *flags,
            ]
        )
        check_file(f"{basename}_Transcript_Quantification.tsv")
    return True


def rewrite_comet_params(source: str | Path, target: str | Path, threads: int) -> None:
    overrides = {
        "decoy_search": 1,
        "output_pepxmlfile": 0,
        "output_percolatorfile": 1,
        "decoy_prefix": "DECOY_",
        "num_threads": threads,
    }
    remaining = dict(overrides)
    lines = []
    for line in Path(source).read_text(encoding="utf-8").splitlines(True):
        key = line.partition("=")[0].strip()
        if key not in remaining:
            lines.append(line)
            continue
        _, _, tail = line.partition("=")
        _, marker, comment = tail.partition("#")
        suffix = f"  #{comment}" if marker else ""
        lines.append(f"{key} = {remaining.pop(key)}{suffix}\n")
    for key, value in remaining.items():
        lines.append(f"{key} = {value}\n")
    Path(target).write_text("".join(lines), encoding="utf-8")


def proteomics(settings: dict, basename: str, query: Path, output_dir: Path) -> bool:
    section = settings["proteomics"]
    if not section["comet-exe"]:
        return False
    if not query.is_file():
        raise PipelineError(f"missing protein query: {query}")

    output_dir.mkdir(parents=True, exist_ok=True)
    comet_output = output_dir / "COMET_Output"
    percolator_raw = output_dir / "PERCOLATOR_Raw"
    percolator_parsed = output_dir / "PERCOLATOR_Parsed"
    for directory in (comet_output, percolator_raw, percolator_parsed):
        directory.mkdir(exist_ok=True)

    spectra_dir = Path(section["comet-spectrometry"]).absolute()
    spectra = sorted(spectra_dir.glob(f"*.{section['comet-ext']}"))
    if not spectra:
        raise PipelineError(f"no Comet input spectra found in {spectra_dir}")

    with contextlib.chdir(output_dir):
        params = Path("comet.params").absolute()
        rewrite_comet_params(
            section["comet-params"], params, settings["pipeline"]["threads"]
        )
        command = [
            str(section["comet-exe"]),
            f"-P{params}",
            f"-D{query}",
        ]
        if section.get("first") is not None:
            command.extend([f"-F{section['first']}", f"-L{section['last']}"])
        command.extend(str(spectrum) for spectrum in spectra)
        run(command)

        pins = {}
        for candidate in (
            *spectra_dir.glob("*.pin"),
            *Path.cwd().glob("*.pin"),
        ):
            pins[candidate.name] = candidate
        if not pins:
            raise PipelineError("Comet produced no PIN files")

        percolator = str(section["percolator-exe"])
        for name, pin in sorted(pins.items()):
            if pin.resolve() != (comet_output / name).resolve():
                shutil.copyfile(pin, comet_output / name)
            stem = re.sub(r"\.[0-9].*", "", pin.stem)
            with contextlib.chdir(percolator_raw):
                run(
                    [
                        percolator,
                        "-r",
                        f"{stem}_peptide_output.tsv",
                        "-m",
                        f"{stem}_percolator_output.tsv",
                        "-B",
                        f"{stem}_decoy_output.tsv",
                        str(pin.absolute()),
                    ]
                )
            with contextlib.chdir(percolator_parsed):
                run(
                    [
                        settings["pipeline"]["python"],
                        str(script("percolator_parser.py")),
                        "-p",
                        str(percolator_raw / f"{stem}_percolator_output.tsv"),
                        "-qv",
                        str(section["percolator-qvalue"]),
                        "-b",
                        f"{basename}_{stem}",
                    ],
                    allowed=(0, 2),
                )

        with contextlib.chdir(percolator_parsed):
            quantitative_proteomics(Path.cwd(), basename)
        sort_tsv(
            percolator_parsed / f"{basename}_pre_total_Proteomics_Quantification.tsv",
            output_dir / f"{basename}_Total_Proteomics_Quantification.tsv",
        )
    return True


def summary(
    settings: dict,
    basename: str,
    rna: bool,
    ms: bool,
    functional_dir: Path,
    transcript_dir: Path,
    peptide_dir: Path,
) -> None:
    command = [
        settings["pipeline"]["python"],
        str(script("summary_parser.py")),
        "-b",
        basename,
        "-annot",
        str(Path(ANNOTATIONS).absolute()),
        "-ipr_annot",
        str(functional_dir / f"{basename}_interproscan_annotated_output.gff3"),
        "-ipr_hyp",
        str(functional_dir / f"{basename}_interproscan_hypothetical_output.gff3"),
    ]
    if rna:
        command.extend(
            [
                "-tr",
                str(transcript_dir / f"{basename}_Transcript_Quantification.tsv"),
            ]
        )
    if ms:
        command.extend(
            [
                "-proteomics",
                str(peptide_dir / f"{basename}_Total_Proteomics_Quantification.tsv"),
            ]
        )
    run(command)


@app.command()
def main(
    seq: str | None = typer.Option(None, "-s", "--seq"),
    protein: str | None = typer.Option(None, "-p", "--prot"),
    config: str = typer.Option(..., "-c", "--config"),
    gff: str | None = typer.Option(None, "-gff"),
) -> None:
    """Run the legacy pipeline from a TOML configuration file."""
    if (seq is None) == (protein is None):
        raise typer.BadParameter("exactly one of -s/--seq or -p/--prot is required")
    source_input = seq or protein
    if source_input is None:  # unreachable after validation; narrows the type
        raise typer.BadParameter("exactly one of -s/--seq or -p/--prot is required")
    if gff and protein is None:
        raise typer.BadParameter("-gff requires -p/--prot")

    configure_logging()
    try:
        settings = load(config, protein_mode=protein is not None)
    except ConfigError as exc:
        LOGGER.error(f"{exc}")
        raise typer.BadParameter(str(exc)) from exc

    absolute_source = Path(source_input).absolute()
    inputs: dict[str, str | None] = {}
    for key, value in (("seq", seq), ("protein", protein), ("gff", gff)):
        inputs[key] = str(Path(value).absolute()) if value else None
    for key, value in inputs.items():
        if value and not Path(value).is_file():
            raise typer.BadParameter(f"{key} input is not a file: {value}")

    basename = settings["pipeline"]["basename"]
    home = Path(f"AnnotaPipeline_{basename}").absolute()
    home.mkdir(exist_ok=True)
    augustus_dir = home / f"1_GenePrediction_{basename}"
    similarity_dir = home / f"2_SimilarityAnalysis_{basename}"
    functional_dir = home / f"3_FunctionalAnnotation_{basename}"
    transcript_dir = home / f"4_TranscriptQuantification_{basename}"
    peptide_index = "5" if seq else "4"
    peptide_dir = home / f"{peptide_index}_PeptideIdentification_{basename}"
    for directory in (augustus_dir, similarity_dir, functional_dir):
        directory.mkdir(exist_ok=True)

    augustus_gff: Path | None = None
    clean_fasta = ""
    transcripts: Path | None = None

    try:
        os.chdir(augustus_dir)
        source = absolute_source
        source_fasta = source.name
        shutil.copyfile(source, source_fasta)
        if seq:
            augustus_bin, augustus_config, augustus_script = augustus_paths(settings)
            augustus_gff = augustus_dir / f"AUGUSTUS_{basename}.gff"
            with augustus_gff.open("w", encoding="utf-8") as output:
                run(
                    [
                        str(augustus_bin),
                        f"--AUGUSTUS_CONFIG_PATH={augustus_config}",
                        *augustus_flags(settings),
                        source_fasta,
                    ],
                    stdout=output,
                )
            check_file(augustus_gff)
            run(
                [
                    str(augustus_script),
                    augustus_gff.name,
                    f"--seqfile={source}",
                ]
            )
            source_fasta = f"AUGUSTUS_{basename}.aa"
            check_file(source_fasta)
            transcripts = augustus_dir / f"AUGUSTUS_{basename}.codingseq"
            check_file(transcripts)
        sequence_cleaner(
            source_fasta,
            settings["seq_cleaner"]["minsize-seq"],
            f"Clear_{basename}.fasta",
        )
        clean_fasta = str((augustus_dir / f"Clear_{basename}.fasta").absolute())

        os.chdir(similarity_dir)
        secondary = settings["databases"]
        db_flag = {
            "eupathdb": "-spdb",
            "nrdb": "-nr",
            "trembldb": "-trbl",
            "custom": "-spdb",
        }[secondary["secondary-format"]]
        blast_opts = settings["local_aligner"]
        run(
            [
                settings["pipeline"]["python"],
                str(script("blastp_parser.py")),
                "-s",
                clean_fasta,
                "-sp",
                secondary["swissprot-db"],
                "-basename",
                basename,
                db_flag,
                secondary["secondary-db"],
                "-id",
                str(blast_opts["identity"]),
                "-pos",
                str(blast_opts["positivity"]),
                "-cov",
                str(blast_opts["coverage"]),
                "-kw",
                ",".join(settings["pipeline"]["keywords"]),
                "-t",
                str(settings["pipeline"]["threads"]),
                "-max_target_seqs",
                str(blast_opts["max_target_seqs"]),
                "-evalue",
                str(blast_opts["evalue"]),
                "-blastp",
                blast_opts["blast-exe"],
                "-aligner",
                blast_opts["aligner"],
                "-customsep",
                secondary["customsep"],
                "-customcolumn",
                str(secondary["customcolumn"]),
            ]
        )
        for suffix in ("annotated", "hypothetical", "no_hit"):
            check_file(f"{basename}_{suffix}_products.txt")

        os.chdir(functional_dir)
        split_blast_products(basename, clean_fasta, similarity_dir)
        for kind in ("hypothetical", "annotated"):
            interproscan(settings, kind, basename)
        hmmer = settings["hmmer"]
        hmm_out = f"{basename}_hmmscan_output.txt"
        run(
            [
                hmmer["hmmscan-exe"],
                "--cpu",
                str(settings["pipeline"]["threads"]),
                "--tblout",
                hmm_out,
                "--noali",
                "-E",
                str(hmmer["E"]),
                "--domE",
                str(hmmer["domE"]),
                *tool_flags(hmmer, skip={"hmmscan-exe", "E", "domE"}),
                secondary["pfam"],
                "Hypothetical_Products.fasta",
            ],
            stdout=subprocess.DEVNULL,
        )
        check_file(hmm_out, empty=True)
        rps = settings["rpsblast"]
        run(
            [
                rps["rpsblast-exe"],
                "-query",
                "Hypothetical_Products.fasta",
                "-out",
                f"{basename}_rpsblast_output.tsv",
                "-db",
                secondary["cdd-db"],
                "-outfmt",
                BLAST_FIELDS,
                "-num_threads",
                str(settings["pipeline"]["threads"]),
                "-evalue",
                str(rps["evalue"]),
                "-max_target_seqs",
                str(rps["max_target_seqs"]),
                *tool_flags(
                    rps,
                    skip={
                        "rpsblast-exe",
                        "evalue",
                        "max_target_seqs",
                    },
                ),
            ]
        )
        check_file(f"{basename}_rpsblast_output.tsv")
        make_annotations(settings, basename, similarity_dir, functional_dir, home)

        os.chdir(home)
        export_sequences(
            settings,
            basename,
            inputs["protein"],
            inputs["gff"],
            str(augustus_gff) if augustus_gff else None,
            clean_fasta,
        )
        query = home / f"AnnotaPipeline_{basename}_proteins.fasta"
        rna_done = kallisto(settings, basename, transcripts, transcript_dir)
        ms_done = proteomics(settings, basename, query, peptide_dir)
        summary(
            settings,
            basename,
            rna_done,
            ms_done,
            functional_dir,
            transcript_dir,
            peptide_dir,
        )
    except (PipelineError, OSError) as exc:
        LOGGER.error(f"{exc}")
        raise typer.Exit(1) from exc
    LOGGER.info(f"AnnotaPipeline completed: {home}")


if __name__ == "__main__":
    app()
