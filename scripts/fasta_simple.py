#!/usr/bin/env python3
"""Legacy CLI for joining FASTA sequences with annotations."""

import warnings

import typer
from sequence_utils import read_annotation_table

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    annot: str = typer.Option(..., "-annot", "--annotations"),
    basename: str = typer.Option(..., "-b", "--basename"),
    seqs: str = typer.Option(..., "-faf", "--fasta"),
    org: str = typer.Option(..., "-org", "--organism"),
) -> None:
    ant = read_annotation_table(annot)

    from Bio import SeqIO

    # ============== Write fasta file with new header ==========================

    ids_warn = []
    with open(f"AnnotaPipeline_{basename!s}_proteins.fasta", "w") as out:
        for record in SeqIO.parse(seqs, "fasta"):
            if record.id not in ant:
                ids_warn.append(record.id)
                continue
            record.description = f"| Organism: {org} | Description: {ant[record.id]}"
            SeqIO.write([record], out, "fasta-2line")

    # =============== Warning message, if something goes wrong =====================

    if len(ids_warn) > 0:
        warnings.warn(
            "WARNING: Not all sequences from fasta file were in annotation file",
            stacklevel=2,
        )
        with open(f"{basename!s}_ids_with_no_annotations.txt", "w") as warn_seq:
            warn_seq.writelines(str(a) + "\n" for a in ids_warn)
        warnings.warn("INFO: IDs stored in ids_with_no_annotations.txt", stacklevel=2)


if __name__ == "__main__":
    app()
