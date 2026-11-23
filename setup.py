from pathlib import Path

from setuptools import find_packages, setup

requirements = [
    "biopython>=1.85",
    "pandas>=2.2.3",
    "typer>=0.12",
]

long_description = Path("README.md").read_text(encoding="utf-8")

setup(
    name="AnnotaPipeline",
    packages=find_packages(),
    py_modules=[
        "annotation_io",
        "blast_classification",
        "config",
        "sequence_utils",
    ],
    package_dir={"": "scripts"},
    version="2.0.0",
    description="Pipeline to predict and annotate eukaryotic genes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    author="Guilherme Augusto Maia, Eric Kazuo Kawagoe, Vilmar Benetti Filho, "
    "Tatiany Aparecida Soratto, Renato Moreira Simões, Glauber Wagner",
    author_email="labinfo.ufsc@gmail.com",
    license="Apache-2.0",
    zip_safe=False,
    scripts=[
        "scripts/annotapipeline.py",
        "scripts/blastp_parser.py",
        "scripts/fasta_simple.py",
        "scripts/fastatogff.py",
        "scripts/funcannotation_parser.py",
        "scripts/gfftofasta_parser.py",
        "scripts/info_parser.py",
        "scripts/kallisto_parser.py",
        "scripts/percolator_parser.py",
        "scripts/summary_parser.py",
    ],
    include_package_data=True,
    python_requires=">=3.12",
    url="https://github.com/labioinfoufsc/AnnotaPipeline",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.12",
    ],
)
