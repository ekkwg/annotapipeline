#!/usr/bin/env python3
"""Legacy CLI for joining annotations into a GFF file."""

import re

import typer
from sequence_utils import read_annotation_table

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    anot: str = typer.Option(..., "-all", "--allah"),
    gff: str = typer.Option(..., "-gff", "--gff"),
    basename: str = typer.Option(..., "-b", "--basename"),
) -> None:
    def transcript_id(fields):
        value = re.sub(r";.*", "", fields[-1])
        return re.sub(r".*=", "", value)

    def write_block(block, annotations, outfile):
        for line in block:
            identifier = None
            newline = ""
            if "\ttranscript\t" in line:
                newline = "\n" if line.endswith("\n") else ""
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 9:
                    outfile.write(line)
                    continue
                identifier = transcript_id(fields)
                if identifier not in annotations:
                    outfile.write(line)
                    continue
                attributes = fields[8].rstrip(";").split(";")
                annotation = annotations[identifier].replace("\t", " ").strip()
                attributes.append(f"annotation={annotation}")
                fields[8] = ";".join(attributes)
                outfile.write("\t".join(fields) + newline)
            else:
                outfile.write(line)

    dict_anot = read_annotation_table(anot)
    with (
        open(str(gff)) as source,
        open(f"{basename!s}_Annotated_GFF.gff", "w") as outfile,
    ):
        block = []
        started = False
        for line in source:
            if line.startswith("# start gene "):
                if block:
                    write_block(block, dict_anot, outfile)
                block = [line]
                started = True
            elif not started:
                outfile.write(line)
            else:
                block.append(line)
        if block:
            write_block(block, dict_anot, outfile)


if __name__ == "__main__":
    app()
