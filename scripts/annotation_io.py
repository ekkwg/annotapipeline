"""Sequence preparation and proteomics aggregation shared by both runners."""

import logging
from pathlib import Path


def sequence_cleaner(fasta_file, min_length=0, output_file=None):
    """Write records meeting the minimum length, preserving input order."""
    from Bio import SeqIO

    destination = output_file or f"Clear_{Path(fasta_file).name}"
    with open(destination, "w") as output:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if len(record.seq) >= min_length:
                output.write(f">{record.id}\n{str(record.seq).upper()}\n")
    return Path(destination)


def fasta_fetcher(input_fasta, id_list, fetcher_output):
    """Write FASTA records whose identifiers exactly match ``id_list``."""
    from Bio import SeqIO

    wanted = set(id_list)
    with SeqIO.parse(input_fasta, "fasta") as source:
        records = (record for record in source if record.id in wanted)
        written = SeqIO.write(records, fetcher_output, "fasta")
    if written != len(wanted):
        logging.getLogger(__name__).warning(
            f"{written} of {len(wanted)} requested FASTA IDs found in {input_fasta}"
        )
    return written


def annotate_codingseq(aa_fasta, codingseq_fasta, basename):
    """Transfer protein annotations to coding sequences by exact ID."""
    from Bio import SeqIO

    with SeqIO.parse(aa_fasta, "fasta") as source:
        annotations = {record.id: record.description for record in source}
    output = Path(f"AnnotaPipeline_{basename}_transcripts.fasta")
    with (
        output.open("w") as handle,
        SeqIO.parse(codingseq_fasta, "fasta") as source,
    ):
        for record in source:
            if record.id not in annotations:
                continue
            record.description = annotations[record.id]
            record.seq = record.seq.upper()
            SeqIO.write(record, handle, "fasta-2line")
    return output


def add_features(feature, data, save):
    """Compute the legacy unique and total counts for one feature."""
    unique_column = f"Unique {feature}"
    total_column = f"Total {feature}"
    selected = data.copy()
    singletons = selected.drop_duplicates(feature, keep=False)
    selected[unique_column] = selected.index.isin(singletons.index)
    selected.replace({False: 0, True: 1}, inplace=True)
    totals = (
        selected.groupby("ProteinID")
        .size()
        .sort_values(ascending=False)
        .reset_index(name=total_column)
    )
    unique = (
        selected.loc[selected[unique_column] == 1]
        .drop(columns=["Peptide", "Spectrum"])
        .groupby("ProteinID")
        .size()
        .reset_index(name=unique_column)
    )
    if "ProteinID" not in save.columns:
        save["ProteinID"] = totals["ProteinID"]
    save = save.set_index("ProteinID").join(unique.set_index("ProteinID"))
    return save.join(totals.set_index("ProteinID")).reset_index()


def aggregate_proteomics(parsed_files, output_file):
    import pandas as pd

    data = pd.DataFrame({"ProteinID": [], "Peptide": [], "Spectrum": []})
    for path in parsed_files:
        data = pd.concat([data, pd.read_csv(path, sep="\t", header=0)])

    # retain legacy counting, including indices and spectrum identifiers
    total = add_features("Peptide", data, pd.DataFrame())
    total = add_features("Spectrum", data, total)
    total = total.fillna(0).astype({"Unique Peptide": int, "Unique Spectrum": int})
    total.sort_values(by="ProteinID", ascending=False).to_csv(
        output_file, sep="\t", index=False
    )


def quantitative_proteomics(path, basename):
    aggregate_proteomics(
        Path(path).glob("*_parsed.tsv"),
        f"{basename}_pre_total_Proteomics_Quantification.tsv",
    )
