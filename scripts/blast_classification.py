"""BLAST/DIAMOND invocation and annotation classification shared by runners."""

from __future__ import annotations

import re
import subprocess
from itertools import groupby
from pathlib import Path
from typing import Any

from Bio import SeqIO

BLAST_FIELDS = "6 qseqid sseqid sacc bitscore evalue ppos pident qcovs stitle"
DIAMOND_FIELDS = "6 qseqid sseqid sseqid bitscore evalue ppos pident qcovhsp stitle"
ALIGNERS = ("blast", "diamond")


def blast_command(
    executable: str | Path,
    query: str | Path,
    output: str | Path,
    database: str | Path,
    max_targets: int,
    evalue: float,
    cpus: int,
    aligner: str = "blast",
) -> list[str]:
    """Build the tabular search command for BLAST+ or DIAMOND."""
    cpus = max(int(cpus), 1)
    if aligner == "blast":
        return [
            str(executable),
            "-query",
            str(query),
            "-out",
            str(output),
            "-db",
            str(database),
            "-evalue",
            str(evalue),
            "-outfmt",
            BLAST_FIELDS,
            "-max_target_seqs",
            str(max_targets),
            "-num_threads",
            str(cpus),
        ]
    if aligner == "diamond":
        # DIAMOND requires outfmt fields as separate arguments and lacks
        # BLAST's sacc/qcovs fields; sseqid fills the unused accession column
        # and qcovhsp fills the coverage column.
        return [
            str(executable),
            "blastp",
            "--query",
            str(query),
            "--out",
            str(output),
            "--db",
            str(database),
            "--evalue",
            str(evalue),
            "--outfmt",
            *DIAMOND_FIELDS.split(),
            "--max-target-seqs",
            str(max_targets),
            "--threads",
            str(cpus),
        ]
    raise ValueError(f"unsupported aligner: {aligner}")


def blast(
    executable: str | Path,
    query: str | Path,
    output: str | Path,
    database: str | Path,
    max_targets: int,
    evalue: float,
    cpus: int,
    aligner: str = "blast",
) -> None:
    """Run a tabular protein search, writing an empty table for no query."""
    if not Path(query).stat().st_size:
        Path(output).write_text("")
        return
    command = blast_command(
        executable,
        query,
        output,
        database,
        max_targets,
        evalue,
        cpus,
        aligner,
    )
    subprocess.run(command, check=True)
    if not Path(output).is_file():
        raise FileNotFoundError(str(output))


def description(
    title: str, db_format: str, separator: str = "|", column: int = 5
) -> str:
    """Extract a product description from a subject title."""
    if db_format == "nrdb":
        match = re.search(r"\s(.*?)\s\[.*", title)
        return (match.group(1) if match else title).strip()
    if db_format in ("swissprot", "trembldb"):
        return title.split(" ", 1)[-1].split("OS=")[0].strip()
    if db_format in ("eupathdb", "custom"):
        return title.split(separator)[column].replace("transcript_product=", "").strip()
    raise ValueError(f"unsupported database format: {db_format}")


def candidates(hits: str | Path, db_format: str, settings: dict[str, Any]):
    """Yield accepted (bitscore, description) pairs grouped by query ID."""
    keywords = [word.lower() for word in settings["keywords"]]

    def rows(source):
        for number, line in enumerate(source, 1):
            fields = line.rstrip("\r\n").split("\t", 8)
            if len(fields) != 9 or not fields[0]:
                raise ValueError(f"{hits}:{number}: invalid tabular hit row")
            yield fields

    with open(hits) as source:
        for query, group in groupby(rows(source), key=lambda fields: fields[0]):
            accepted = []
            for fields in group:
                title = description(
                    fields[8],
                    db_format,
                    settings.get("customsep", "|"),
                    settings.get("customcolumn", 5),
                )
                if (
                    float(fields[7]) > float(settings["coverage"])
                    and not any(word in title.lower() for word in keywords)
                    and float(fields[5]) >= float(settings["positivity"])
                    and float(fields[6]) >= float(settings["identity"])
                ):
                    accepted.append((float(fields[3]), title))
            yield query, accepted


def write_candidate(query, accepted, products, report) -> None:
    # max retains the first hit for equal bit scores
    best = max(accepted, key=lambda hit: hit[0])[1]
    titles = ";".join(title for _, title in accepted)
    products.write(f"{query}\t{best}\n")
    report.write(f"{query}\t{len(accepted)} Annotation(s): [{titles}]\n")


def process_swiss(
    basename: str, protein_seq: str | Path, hits: str | Path, settings: dict
) -> None:
    """Select accepted SwissProt products and write the unresolved FASTA."""
    accepted_ids = []
    with (
        open(f"{basename}_annotated_products.txt", "w") as products,
        open(f"{basename}_SwissProt_annotations.txt", "w") as report,
    ):
        for query, accepted in candidates(hits, "swissprot", settings):
            if accepted:
                write_candidate(query, accepted, products, report)
                accepted_ids.append(query)

    accepted = set(accepted_ids)
    unresolved = [
        record
        for record in SeqIO.parse(protein_seq, "fasta")
        if record.id not in accepted
    ]
    with open(f"{basename}_BLASTp_AA_SwissProted.fasta", "w") as output:
        SeqIO.write(unresolved, output, "fasta")


def process_secondary(
    basename: str,
    unresolved: str | Path,
    hits: str | Path,
    db_format: str,
    settings: dict,
) -> None:
    """Classify second-pass hits and write the no-hit protein list."""
    seen = set()
    with (
        open(f"{basename}_hypothetical_products.txt", "w") as hypothetical,
        open(f"{basename}_annotated_products.txt", "a") as products,
        open(f"{basename}_SpecifiedDB_annotations.txt", "w") as report,
    ):
        for query, accepted in candidates(hits, db_format, settings):
            seen.add(query)
            if accepted:
                write_candidate(query, accepted, products, report)
            else:
                hypothetical.write(f"{query}\n")

    with (
        open(unresolved) as fasta,
        open(f"{basename}_no_hit_products.txt", "w") as no_hits,
    ):
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id not in seen:
                no_hits.write(record.id + "\n")
