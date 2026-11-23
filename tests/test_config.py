"""Configuration loader contract checks."""

import os
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from config import ConfigError, as_flag, load  # noqa: E402

MINIMAL_TOML = """
[pipeline]
basename = "Test"
organism = "Model plant"

[databases]
swissprot-db = "swiss.fasta"
secondary-db = "secondary.fasta"
secondary-format = "nrdb"

[local_aligner]
blast-exe = "blastp"
"""


class ConfigTest(unittest.TestCase):
    def write(self, root: Path, name: str, content: str) -> Path:
        path = root / name
        path.write_text(content, encoding="utf-8")
        return path

    def test_toml_defaults_and_explicit_types(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp), "config.toml", MINIMAL_TOML)
            config = load(path, protein_mode=True)
        self.assertEqual(config["pipeline"]["threads"], 1)
        self.assertEqual(config["pipeline"]["keywords"][0], "fragment")
        self.assertEqual(config["local_aligner"]["aligner"], "blast")
        self.assertEqual(config["local_aligner"]["identity"], 40)
        self.assertEqual(config["databases"]["customcolumn"], 5)
        self.assertFalse(config["kallisto"]["single_end"])

    def test_diamond_aligner_is_accepted(self):
        content = MINIMAL_TOML.replace(
            '[local_aligner]\nblast-exe = "blastp"',
            '[local_aligner]\naligner = "diamond"\nblast-exe = "diamond"',
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp), "config.toml", content)
            config = load(path, protein_mode=True)
        self.assertEqual(config["local_aligner"]["aligner"], "diamond")

    def test_unknown_aligner_is_rejected(self):
        content = MINIMAL_TOML.replace(
            '[local_aligner]\nblast-exe = "blastp"',
            '[local_aligner]\naligner = "mmseqs"\nblast-exe = "mmseqs"',
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp), "config.toml", content)
            with self.assertRaisesRegex(ConfigError, "aligner"):
                load(path, protein_mode=True)

    def test_non_toml_suffix_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp), "config.yaml", MINIMAL_TOML)
            with self.assertRaisesRegex(ConfigError, "must use the .toml"):
                load(path, protein_mode=True)

    def test_scheduler_threads_from_environment(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp), "config.toml", MINIMAL_TOML)
            with unittest.mock.patch.dict(
                os.environ, {"SLURM_CPUS_PER_TASK": "7"}, clear=False
            ):
                config = load(path, protein_mode=True)
        self.assertEqual(config["pipeline"]["threads"], 7)

    def test_explicit_threads_override_scheduler(self):
        content = MINIMAL_TOML.replace(
            'organism = "Model plant"',
            'organism = "Model plant"\nthreads = 3',
        )
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp), "config.toml", content)
            with unittest.mock.patch.dict(
                os.environ, {"SLURM_CPUS_PER_TASK": "7"}, clear=False
            ):
                config = load(path, protein_mode=True)
        self.assertEqual(config["pipeline"]["threads"], 3)

    def test_genome_mode_requires_augustus_settings(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp), "config.toml", MINIMAL_TOML)
            with self.assertRaisesRegex(ConfigError, "augustus.strand"):
                load(path, protein_mode=False)

    def test_invalid_basename_is_rejected(self):
        content = MINIMAL_TOML.replace('basename = "Test"', 'basename = "a/b"')
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(Path(tmp), "config.toml", content)
            with self.assertRaisesRegex(ConfigError, "basename"):
                load(path, protein_mode=True)

    def test_unknown_section_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.write(
                Path(tmp), "config.toml", MINIMAL_TOML + "\n[typo]\na = 1\n"
            )
            with self.assertRaisesRegex(ConfigError, "unknown configuration"):
                load(path, protein_mode=True)

    def test_missing_file_is_reported(self):
        with self.assertRaisesRegex(ConfigError, "cannot read configuration"):
            load("does-not-exist.toml", protein_mode=True)

    def test_flag_values(self):
        self.assertEqual(as_flag("on"), "on")
        self.assertEqual(as_flag("Flag"), "flag")
        with self.assertRaises(ConfigError):
            as_flag("yes")

    def test_shipped_example_loads(self):
        load(ROOT / "config" / "config_example.toml", protein_mode=True)

    def test_template_requires_placeholder_replacement(self):
        with self.assertRaisesRegex(ConfigError, "pipeline.basename is required"):
            load(ROOT / "config" / "annotapipeline.toml", protein_mode=True)


if __name__ == "__main__":
    unittest.main()
