"""Canonical schema, defaults, validation and loading for the legacy runner."""

from __future__ import annotations

import os
import re
import tomllib
from pathlib import Path
from typing import Any

DEFAULT_KEYWORDS = (
    "fragment",
    "hypothetical",
    "partial",
    "uncharacterized",
    "unknown",
    "unspecified",
)

_SECTIONS = {
    "pipeline",
    "databases",
    "augustus",
    "augustus_optional",
    "seq_cleaner",
    "local_aligner",
    "rpsblast",
    "interproscan",
    "hmmer",
    "kallisto",
    "proteomics",
}
_REQUIRED = {
    "pipeline": ("basename", "organism"),
    "databases": ("swissprot-db", "secondary-db", "secondary-format"),
    "local_aligner": ("blast-exe",),
}
_ALIGNERS = {"blast", "diamond"}
_SECONDARY_FORMATS = {"eupathdb", "nrdb", "trembldb", "custom"}
_ON_OFF = {"on", "off"}
_BASENAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
_NUMERIC = {
    "threads",
    "customcolumn",
    "minsize-seq",
    "identity",
    "positivity",
    "coverage",
    "max_target_seqs",
    "E",
    "domE",
    "bootstrap",
    "l",
    "s",
    "first",
    "last",
}
_BOOLEAN_TEXT = {"flag", "on", "off"}


class ConfigError(ValueError):
    """Raised when a configuration file is invalid."""


def _default() -> dict[str, Any]:
    return {
        "pipeline": {
            "python": "python3",
            "perl": "perl",
            "threads": _scheduler_threads(),
            "organism": None,
            "basename": None,
            "keywords": list(DEFAULT_KEYWORDS),
        },
        "databases": {
            "swissprot-db": None,
            "secondary-db": None,
            "secondary-format": None,
            "customsep": "|",
            "customcolumn": 5,
            "cdd-db": None,
            "pfam": None,
        },
        "augustus": {
            "augustus-path": "conda",
            "strand": None,
            "genemodel": None,
            "species": None,
            "protein": None,
            "introns": None,
            "start": None,
            "stop": None,
            "cds": None,
        },
        "augustus_optional": {
            "hintsfile": None,
            "extrinsicCfgFile": None,
            "UTR": None,
        },
        "seq_cleaner": {"minsize-seq": 0},
        "local_aligner": {
            "aligner": "blast",
            "blast-exe": "blastp",
            "evalue": 1e-5,
            "identity": 40,
            "positivity": 60,
            "coverage": 30,
            "max_target_seqs": 10,
        },
        "rpsblast": {
            "rpsblast-exe": "rpsblast",
            "evalue": 1e-5,
            "max_target_seqs": 10,
        },
        "interproscan": {"interproscan-exe": "interproscan.sh", "ms": None},
        "hmmer": {
            "hmmscan-exe": "hmmscan",
            "E": 1e-5,
            "domE": 1e-5,
        },
        "kallisto": {
            "kallisto-exe": None,
            "rna-seq": None,
            "l": None,
            "s": None,
            "threshold": None,
            "bootstrap": 0,
        },
        "proteomics": {
            "comet-exe": None,
            "comet-params": None,
            "comet-spectrometry": None,
            "comet-ext": None,
            "first": None,
            "last": None,
            "percolator-exe": "percolator",
            "percolator-qvalue": 0.05,
        },
    }


def _scheduler_threads() -> int:
    for name in ("SLURM_CPUS_PER_TASK", "NSLOTS"):
        raw = os.environ.get(name)
        if raw and raw.isdigit() and int(raw) > 0:
            return int(raw)
    return 1


def _parse(path: Path) -> dict[str, Any]:
    if path.suffix.lower() != ".toml":
        raise ConfigError(f"{path}: configuration must use the .toml suffix")
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ConfigError(f"cannot read configuration {path}: {exc}") from exc
    try:
        document = tomllib.loads(text)
    except tomllib.TOMLDecodeError as exc:
        raise ConfigError(f"invalid TOML in {path}: {exc}") from exc
    if not isinstance(document, dict):
        raise ConfigError(f"{path} must contain a configuration table")
    return document


def _validate(config: dict[str, Any], path: Path, protein_mode: bool) -> None:
    def fail(message: str) -> None:
        raise ConfigError(f"{path}: {message}")

    pipeline = config["pipeline"]
    for field in _REQUIRED["pipeline"]:
        if not pipeline[field]:
            fail(f"pipeline.{field} is required")
    if not _BASENAME.fullmatch(str(pipeline["basename"])):
        fail(
            "pipeline.basename must start with a letter/digit and use only "
            "letters, digits, '.', '_' or '-'"
        )
    if not isinstance(pipeline["keywords"], list) or any(
        not isinstance(word, str) or not word for word in pipeline["keywords"]
    ):
        fail("pipeline.keywords must be a list of nonempty strings")
    for field in (
        "threads",
        "customcolumn",
        "minsize-seq",
        "max_target_seqs",
        "ms",
        "bootstrap",
        "first",
        "last",
    ):
        section = (
            config["databases"]
            if field == "customcolumn"
            else config["seq_cleaner"]
            if field == "minsize-seq"
            else config["local_aligner"]
            if field == "max_target_seqs"
            else config["interproscan"]
            if field == "ms"
            else config["kallisto"]
            if field == "bootstrap"
            else config["proteomics"]
        )
        value = section.get(field)
        if value is not None and (
            isinstance(value, bool) or not isinstance(value, int) or value < 0
        ):
            fail(f"{field} must be a non-negative integer")
    for section_name, field in (
        ("local_aligner", "identity"),
        ("local_aligner", "positivity"),
        ("local_aligner", "coverage"),
        ("proteomics", "percolator-qvalue"),
    ):
        value = config[section_name].get(field)
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not 0 <= value <= 100
        ):
            fail(f"{field} must be a number between 0 and 100")
    qvalue = config["proteomics"]["percolator-qvalue"]
    if qvalue is not None and qvalue > 1:
        fail("proteomics.percolator-qvalue must be between 0 and 1")

    db = config["databases"]
    for field in _REQUIRED["databases"]:
        if not db[field]:
            fail(f"databases.{field} is required")
    if db["secondary-format"] not in _SECONDARY_FORMATS:
        fail("databases.secondary-format must be eupathdb, nrdb, trembldb or custom")
    if db["secondary-format"] == "custom" and not (
        db["customsep"] and db["customcolumn"] is not None
    ):
        fail("custom databases require customsep and customcolumn")
    if not config["local_aligner"]["blast-exe"]:
        fail("local_aligner.blast-exe is required")
    if config["local_aligner"]["aligner"] not in _ALIGNERS:
        fail(f"local_aligner.aligner must be one of: {', '.join(_ALIGNERS)}")

    for field, value in config["augustus"].items():
        if field == "augustus-path" or value is None or value == "":
            continue
        if field in {"strand", "genemodel"}:
            allowed = (
                {"both", "forward", "backward"}
                if field == "strand"
                else {
                    "partial",
                    "intronless",
                    "complete",
                    "atleastone",
                    "exactlyone",
                }
            )
            if str(value).lower() not in allowed:
                fail(f"augustus.{field} must be one of: {', '.join(sorted(allowed))}")
        elif (
            field in {"protein", "introns", "start", "stop", "cds"}
            and str(value).lower() not in _ON_OFF
        ):
            fail(f"augustus.{field} must be on or off")
    if not protein_mode:
        for field, value in config["augustus"].items():
            if field != "augustus-path" and (value is None or value == ""):
                fail(f"augustus.{field} is required for genome input")
        if not config["augustus"]["species"]:
            fail("augustus.species is required for genome input")

    optional = config["augustus_optional"]
    for field in ("hintsfile", "extrinsicCfgFile"):
        if optional.get(field) and not Path(optional[field]).is_file():
            fail(f"augustus_optional.{field} is not a regular file")
    if optional.get("UTR") and str(optional["UTR"]).lower() not in (_ON_OFF | {"flag"}):
        fail("augustus_optional.UTR must be on, off or flag")

    kallisto = config["kallisto"]
    reads = kallisto.get("rna-seq") or []
    if not isinstance(reads, list) or any(not str(item) for item in reads):
        fail("kallisto.rna-seq must be a list of paths")
    threshold = kallisto["threshold"]
    if threshold is not None and not (
        str(threshold).lower() in {"mean", "median"}
        or (isinstance(threshold, (int, float)) and not isinstance(threshold, bool))
    ):
        fail("kallisto.threshold must be mean, median or a number")
    if kallisto["kallisto-exe"]:
        if not 1 <= len(reads) <= 2:
            fail("kallisto.rna-seq must contain one or two read files")
        if threshold is None:
            fail("kallisto.threshold is required")
        if len(reads) == 1 and (kallisto.get("l") is None or kallisto.get("s") is None):
            fail("single-end Kallisto requires l and s")
    kallisto["single_end"] = len(reads) == 1

    proteomics = config["proteomics"]
    if proteomics["comet-exe"]:
        for field in ("comet-params", "comet-spectrometry", "comet-ext"):
            if not proteomics.get(field):
                fail(f"proteomics.{field} is required with Comet")
        if (proteomics.get("first") is None) != (proteomics.get("last") is None):
            fail("proteomics.first and last must be given together")
        if (
            proteomics.get("first") is not None
            and proteomics["last"] <= proteomics["first"]
        ):
            fail("proteomics.last must be greater than first")
        if not proteomics.get("percolator-exe"):
            fail("proteomics.percolator-exe is required with Comet")


def load(path: str | Path, *, protein_mode: bool = False) -> dict[str, Any]:
    """Load and validate a TOML legacy-runner configuration."""
    source = Path(path)
    document = _parse(source)
    sections = {}
    for name, section in document.items():
        if name not in _SECTIONS:
            raise ConfigError(f"{source}: unknown configuration section: {name}")
        if not isinstance(section, dict):
            raise ConfigError(f"{source}: section {name} must be a table")
        sections[name] = section
    config = _default()
    for name, values in sections.items():
        allowed = config[name]
        for field, value in values.items():
            if field not in allowed:
                raise ConfigError(f"{source}: unknown setting: {name}.{field}")
            allowed[field] = value

    # Zero is a documented TOML sentinel for scheduler-derived CPU counts.
    threads = config["pipeline"]["threads"]
    if threads == 0:
        config["pipeline"]["threads"] = _scheduler_threads()

    _validate(config, source, protein_mode)
    config["config_path"] = str(source)
    return config


def as_flag(value: Any) -> str:
    """Normalise legacy on/off and flag values for command-line use."""
    text = str(value).lower()
    if text not in _BOOLEAN_TEXT:
        raise ConfigError(f"unsupported legacy flag value: {value!r}")
    return text
