"""Fixture-backed end-to-end checks for the legacy runner."""

import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(ROOT / "tests"))

from tool_fixture import PROTEINS  # noqa: E402


class RunnerTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="annota-runner-")
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.protein = self.root / "plant protein.fasta"
        self.protein.write_text(PROTEINS)
        self.tools = self.root / "tools"
        self.tools.mkdir()
        self.config = self.root / "annotapipeline.toml"
        self.settings: dict[str, Any] = {
            "pipeline": {
                "python": sys.executable,
                "threads": 1,
                "organism": "Model plant",
                "basename": "fixture",
                "keywords": ["hypothetical", "fragment", "partial"],
            },
            "databases": {
                "swissprot-db": self.root / "swiss",
                "secondary-db": self.root / "secondary",
                "secondary-format": "nrdb",
                "customsep": "|",
                "customcolumn": 5,
                "pfam": self.root / "pfam.hmm",
                "cdd-db": self.root / "cdd",
            },
            "augustus": {
                "augustus-path": "conda",
                "strand": "both",
                "genemodel": "partial",
                "species": "model_plant",
                "protein": "on",
                "introns": "on",
                "start": "on",
                "stop": "on",
                "cds": "on",
            },
            "augustus_optional": {},
            "seq_cleaner": {"minsize-seq": 0},
            "local_aligner": {
                "aligner": "blast",
                "blast-exe": self.wrapper("blastp"),
            },
            "rpsblast": {"rpsblast-exe": self.wrapper("rpsblast")},
            "interproscan": {"interproscan-exe": self.wrapper("interproscan")},
            "hmmer": {"hmmscan-exe": self.wrapper("hmmscan")},
            "kallisto": {"kallisto-exe": "", "rna-seq": []},
            "proteomics": {"comet-exe": "", "percolator-exe": ""},
        }
        for name in ("swiss", "secondary", "pfam.hmm", "cdd"):
            (self.root / name).write_text("fixture reference\n")
        self.out = self.root / "AnnotaPipeline_fixture"
        self.log = []

    def wrapper(self, tool, output=None):
        output = output or self.tools / tool
        command = [
            sys.executable,
            str(ROOT / "tests" / "tool_fixture.py"),
            tool,
        ]
        output.write_text(f'#!/bin/bash\nexec {" ".join(map(str, command))} "$@"\n')
        output.chmod(0o755)
        return output

    def write_config(self):
        lines = []
        for section, values in self.settings.items():
            lines.append(f"[{section}]")
            for key, value in values.items():
                if isinstance(value, Path):
                    value = str(value)
                lines.append(f"{key} = {json.dumps(value)}")
            lines.append("")
        self.config.write_text("\n".join(lines), encoding="utf-8")

    def invoke(self, *args):
        self.write_config()
        env = dict(os.environ, CONDA_PREFIX=str(self.root))
        env["PATH"] = f"{self.tools}{os.pathsep}{env['PATH']}"
        proc = subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "annotapipeline.py"),
                "-c",
                str(self.config),
                *args,
            ],
            cwd=self.root,
            env=env,
            capture_output=True,
            text=True,
            timeout=180,
        )
        report = proc.stdout + proc.stderr
        log = self.root / "AnnotaPipeline_LogFile.log"
        if log.is_file():
            report += "\n" + log.read_text(encoding="utf-8")
        self.log.append(report)
        self.assertEqual(proc.returncode, 0, report)
        return proc

    def annotations(self):
        return (self.out / "All_Annotated_Products.txt").read_text()

    def test_protein_only(self):
        self.invoke("-p", self.protein.name)
        annotations = self.annotations()
        self.assertIn("geneA\tFirst enzyme", annotations)
        self.assertIn("geneB\tCatalytic B", annotations)
        self.assertIn("geneC\thypothetical protein", annotations)
        self.assertIn("geneD\thypothetical protein", annotations)
        proteins = self.out / "AnnotaPipeline_fixture_proteins.fasta"
        self.assertIn("Organism: Model plant", proteins.read_text())
        similarity = self.out / "2_SimilarityAnalysis_fixture"
        self.assertTrue((similarity / "fixture_BLASTp_AAvsSwissProt.tsv").is_file())
        self.assertFalse(list(similarity.glob("*.outfmt6")))
        self.assertTrue((self.out / "1_GenePrediction_fixture").is_dir())
        self.assertTrue((self.out / "2_SimilarityAnalysis_fixture").is_dir())
        self.assertTrue((self.out / "3_FunctionalAnnotation_fixture").is_dir())
        self.assertFalse(list(self.out.glob("*Annotated_GFF.gff")))
        summary = self.out / "AnnotaPipeline_fixture_Summary.tsv"
        self.assertTrue(summary.is_file())

    def test_diamond_aligner(self):
        self.settings["local_aligner"]["aligner"] = "diamond"
        self.settings["local_aligner"]["blast-exe"] = self.wrapper("diamond")
        self.invoke("-p", self.protein.name)
        self.assertIn("geneA\tFirst enzyme", self.annotations())
        similarity = self.out / "2_SimilarityAnalysis_fixture"
        self.assertTrue((similarity / "fixture_BLASTp_AAvsSwissProt.tsv").is_file())
        self.assertFalse(list(similarity.glob("*.outfmt6")))

    def test_protein_with_gff(self):
        gff = subprocess.check_output(
            [self.wrapper("augustus")], cwd=self.root, text=True
        )
        gff_path = self.root / "supplied.gff"
        gff_path.write_text(gff)
        self.invoke("-p", self.protein.name, "-gff", gff_path.name)
        self.assertTrue((self.out / "fixture_Annotated_GFF.gff").is_file())
        self.assertIn("geneA\tFirst enzyme", self.annotations())

    def test_genome_and_optional_omics(self):
        genome = self.root / "genome.fasta"
        genome.write_text(">scaffold\n" + "ATG" * 32 + "\n")
        augustus_root = self.root / "augustus"
        (augustus_root / "bin").mkdir(parents=True)
        (augustus_root / "scripts").mkdir()
        (augustus_root / "config").mkdir()
        self.wrapper("augustus", augustus_root / "bin" / "augustus")
        self.wrapper("getannofasta", augustus_root / "scripts" / "getAnnoFasta.pl")
        self.settings["augustus"]["augustus-path"] = str(augustus_root)
        reads = []
        for name in ("z-R1.fq", "a-R2.fq"):
            read = self.root / name
            read.write_text("@record\nAC\n+\nII\n")
            reads.append(str(read))
        self.settings["kallisto"] = {
            "kallisto-exe": self.wrapper("kallisto"),
            "rna-seq": reads,
            "threshold": "mean",
            "bootstrap": 0,
        }
        spectra_dir = self.root / "spectra"
        spectra_dir.mkdir()
        for name in ("sample.mgf", "no_hits.mgf"):
            (spectra_dir / name).write_text("fixture spectra\n")
        params = self.root / "input comet.params"
        original = "decoy_search = 0\noutput_percolatorfile = 0\n"
        params.write_text(original)
        self.settings["proteomics"] = {
            "comet-exe": self.wrapper("comet"),
            "comet-params": params,
            "comet-spectrometry": spectra_dir,
            "comet-ext": "mgf",
            "percolator-exe": self.wrapper("percolator"),
            "percolator-qvalue": 0.05,
        }

        self.invoke("-s", genome.name)

        self.assertTrue((self.out / "fixture_Annotated_GFF.gff").is_file())
        transcripts = self.out / "AnnotaPipeline_fixture_transcripts.fasta"
        self.assertIn("ATGGCTGCTGCT", transcripts.read_text())
        self.assertTrue((self.out / "4_TranscriptQuantification_fixture").is_dir())
        self.assertTrue((self.out / "5_PeptideIdentification_fixture").is_dir())
        self.assertEqual(params.read_text(), original)
        self.assertFalse(list(self.root.glob("*.pin")))
        summary = self.out / "AnnotaPipeline_fixture_Summary.tsv"
        with summary.open() as handle:
            rows = {
                row["ProteinID"]: row for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(rows["geneA"]["TPM"], "9.0")
        self.assertEqual(rows["geneA"]["Total_Peptide"], "2")


if __name__ == "__main__":
    unittest.main()
