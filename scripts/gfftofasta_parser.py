#!/usr/bin/env python3
"""Legacy CLI for writing GFF information into FASTA headers."""

import warnings

import typer
from sequence_utils import parse_annotations, parse_gff_transcripts

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    gff: str = typer.Option(..., "-gff"),
    annot: str = typer.Option(..., "-annot", "--annotations"),
    basename: str = typer.Option(..., "-b", "--basename"),
    seqs: str = typer.Option(..., "-faf", "--fasta"),
    org: str = typer.Option(..., "-org", "--organism"),
) -> None:
    def fasta_records(path):
        header = None
        lines = []
        with open(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    if header is not None:
                        yield header, lines
                    header = line[1:].rstrip("\r\n")
                    lines = []
                else:
                    lines.append(line)
        if header is not None:
            yield header, lines

    trans = parse_gff_transcripts(str(gff))
    ant = parse_annotations(str(annot))

    # ============== Write fasta file with new header ==========================

    # Remove quotes from organism
    organism = str(org).replace('"', "").replace("'", "")

    with open(f"AnnotaPipeline_{basename!s}_proteins.fasta", "w") as annotafile:
        ids_warn = []  # store ids with no annotation
        for header, sequence in fasta_records(str(seqs)):
            try:
                annotafile.write(
                    f">{header} | "
                    f"Organism: {organism} | "
                    f"Location: {trans[header]['scaff']!s} | "
                    f"Start: {trans[header]['start']!s} | "
                    f"End: {trans[header]['end']!s} | "
                    f"Strand: {trans[header]['strand']!s} | "
                    f"Description: {ant[header]!s}\n"
                )
                annotafile.write("".join(sequence))
            except Exception:
                # warn ignored, user will be informed all at once
                ids_warn.append(header)

    # =============== Warning message, if something goes wrong =====================

    if len(ids_warn) > 0:
        warnings.warn(
            "WARNING: Not all sequences from fasta file were in annotation file",
            stacklevel=2,
        )
        with open(f"{basename!s}_ids_with_no_annotations.txt", "w") as warn_seq:
            for header in ids_warn:
                warn_seq.write(str(header) + "\n")
        warnings.warn("INFO: IDs stored in ids_with_no_annotations.txt", stacklevel=2)


if __name__ == "__main__":
    app()
