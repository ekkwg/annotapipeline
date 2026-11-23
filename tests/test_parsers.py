"""Focused checks for scientific parser invariants."""

import contextlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from annotation_io import annotate_codingseq, fasta_fetcher  # noqa: E402
from blast_classification import (  # noqa: E402
    BLAST_FIELDS,
    blast_command,
)
from sequence_utils import read_annotation_table  # noqa: E402

FASTA = ">tx1 first\nAAAA\n>tx10 second\nCCCC\n"
GFF = """##gff-version 3
# start gene gene1
scaffold\tAUGUSTUS\tgene\t1\t20\t.\t+\t.\tID=gene1
scaffold\tAUGUSTUS\ttranscript\t1\t20\t.\t+\t.\tID=tx1;Parent=gene1
scaffold\tAUGUSTUS\tCDS\t1\t20\t.\t+\t0\tID=tx1.cds;Parent=tx1
# end gene gene1
"""
ANNOTATIONS = "tx1\thypothetical protein\twith tab\n"


class ParserTest(unittest.TestCase):
    def test_fasta_fetcher_uses_exact_identifiers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.fasta"
            target = root / "selected.fasta"
            source.write_text(FASTA, encoding="utf-8")
            written = fasta_fetcher(source, ["tx1"], target)
            self.assertEqual(written, 1)
            self.assertNotIn("tx10", target.read_text(encoding="utf-8"))

    def test_coding_sequences_match_proteins_by_exact_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            proteins = root / "proteins.fasta"
            coding = root / "coding.fasta"
            proteins.write_text(">tx1 protein description\nKKKK\n")
            coding.write_text(">tx1\nAAAA\n>tx10\nCCCC\n")
            with contextlib.chdir(root):
                output = annotate_codingseq(proteins, coding, "test")
                text = output.read_text(encoding="utf-8")
            self.assertIn("protein description", text)
            self.assertNotIn(">tx10", text)

    def test_annotation_table_preserves_tabs_in_descriptions(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "annotations.tsv"
            path.write_text(ANNOTATIONS, encoding="utf-8")
            self.assertEqual(
                read_annotation_table(path),
                {"tx1": "hypothetical protein\twith tab"},
            )

    def test_blast_command_keeps_legacy_tabular_contract(self):
        command = blast_command("blastp", "q.fa", "out.tsv", "db", 10, 1e-5, 4)
        self.assertEqual(command[0], "blastp")
        self.assertEqual(command[command.index("-outfmt") + 1], BLAST_FIELDS)
        self.assertEqual(command[command.index("-num_threads") + 1], "4")
        self.assertNotIn("blastp", command[1:])

    def test_diamond_command_uses_blastp_and_qcovhsp(self):
        command = blast_command(
            "diamond", "q.fa", "out.tsv", "db", 10, 1e-5, 4, "diamond"
        )
        self.assertEqual(command[:2], ["diamond", "blastp"])
        self.assertIn("--query", command)
        self.assertIn("--out", command)
        self.assertIn("qcovhsp", command)
        self.assertEqual(command[command.index("--threads") + 1], "4")
        self.assertEqual(command[command.index("--max-target-seqs") + 1], "10")

    def test_unsupported_aligner_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unsupported aligner"):
            blast_command("tool", "q.fa", "out.tsv", "db", 10, 1e-5, 1, "mmseqs")

    def test_gff_annotation_preserves_transcript_attributes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            gff = root / "input.gff"
            annotations = root / "annotations.tsv"
            gff.write_text(GFF, encoding="utf-8")
            annotations.write_text(ANNOTATIONS, encoding="utf-8")
            proc = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "fastatogff.py"),
                    "-gff",
                    str(gff),
                    "-all",
                    str(annotations),
                    "-b",
                    "test",
                ],
                cwd=root,
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            output = (root / "test_Annotated_GFF.gff").read_text(encoding="utf-8")
            self.assertIn("ID=tx1;Parent=gene1", output)
            self.assertIn("annotation=hypothetical protein with tab", output)
            self.assertIn("ID=tx1.cds;Parent=tx1", output)


if __name__ == "__main__":
    unittest.main()
