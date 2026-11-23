"""Shared GFF3 and annotation-table helpers for the legacy parsers."""

from __future__ import annotations

import re
from collections.abc import Iterator
from pathlib import Path


def read_annotation_table(path: str | Path) -> dict[str, str]:
    """Read a two-column table without truncating descriptions at tabs."""
    result = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        fields = line.split("\t", 1)
        if len(fields) == 2:
            result[fields[0]] = fields[1]
    return result


def parse_annotations(path: str | Path) -> dict[str, str]:
    """Alias for the legacy parser's annotation table helper name."""
    return read_annotation_table(path)


def parse_gff_transcripts(path: str | Path) -> dict[str, dict[str, str]]:
    """Return transcript coordinates from AUGUSTUS-style GFF3."""
    records = {}
    for fields in gff3_records(path):
        if fields[2] != "transcript":
            continue
        identifier = re.sub(r";.*", "", fields[-1])
        identifier = re.sub(r".*=", "", identifier)
        records[identifier] = {
            "scaff": fields[0],
            "start": fields[3],
            "end": fields[4],
            "strand": fields[6],
        }
    return records


def gff3_records(path: str | Path) -> Iterator[list[str]]:
    """Yield non-comment GFF3 records with their nine standard fields."""
    for number, line in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), 1
    ):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"{path}:{number}: expected nine GFF3 fields")
        yield fields


def iter_interpro_groups(path: str | Path) -> Iterator[tuple[str, list[str]]]:
    """Group InterProScan GFF3 annotations by transcript identifier."""
    current: str | None = None
    lines: list[str] = []
    for fields in gff3_records(path):
        identifier = re.sub(r";.*", "", fields[-1])
        identifier = re.sub(r".*=", "", identifier)
        if current is not None and identifier != current:
            yield current, lines
            lines = []
        current = identifier
        lines.append("\t".join(fields))
    if current is not None:
        yield current, lines
