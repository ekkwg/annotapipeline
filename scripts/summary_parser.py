#!/usr/bin/env python3
"""build the AnnotaPipeline per-protein summary table"""

import re

import typer

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    basename: str = typer.Option(..., "-b"),
    annot: str = typer.Option(..., "-annot"),
    ipr_annot: str = typer.Option(..., "-ipr_annot"),
    ipr_hyp: str = typer.Option(..., "-ipr_hyp"),
    tr: str = typer.Option(None, "-tr"),
    proteomics: str = typer.Option(None, "-proteomics"),
) -> None:
    superfamily_dict = {}
    ipr_dict = {}
    go_dict = {}

    def get_interpro_info(arq_entrada):
        with open(arq_entrada) as source:
            entrada = source.read().split("##FASTA")
        interp = entrada[0].split(
            "##sequence-region"
        )  # keep only the annotation section before FASTA
        del interp[0]

        unwanted_db = ["coils", "gene3d", "mobidblite"]

        def create_or_reset_lists():
            ontologia = []
            interpro = []
            superfamily = []
            return ontologia, interpro, superfamily

        def add_to_dicts(id, go, ipr, superfamily):
            # store unique values in the per-protein dictionaries
            superfamily_dict[id] = ",".join(sorted(set(superfamily)))
            ipr_dict[id] = ",".join(sorted(set(ipr)))
            go_dict[id] = ",".join(sorted(go))

        def add_none_to_empty():
            if not ontologia:
                ontologia.append("None")
            if not interpro:
                interpro.append("None")
            if not superfamily:
                superfamily.append("None")

        # collect ontology and domain evidence for each query
        for seq_reg in interp:
            ontologia, interpro, superfamily = create_or_reset_lists()
            seq_reg = seq_reg.splitlines()
            nome_subject = seq_reg[0].split()[0]
            for linha in seq_reg:
                if linha == seq_reg[0]:
                    linha = linha.split(" ")
                    del linha[0]
                elif linha == seq_reg[1]:
                    if linha == seq_reg[-1]:
                        add_none_to_empty()
                        add_to_dicts(nome_subject, ontologia, interpro, superfamily)
                        ontologia, interpro, superfamily = create_or_reset_lists()
                    # skip the sequence-region header
                else:
                    linha = linha.split("\t")
                    nome_subject = linha[0]
                    db = linha[1]
                    if not any(s in db.lower() for s in unwanted_db):
                        anotation = linha[-1].split(";")
                        for field in anotation:
                            if "Ontology" in field:
                                go = (
                                    field.replace('"', "")
                                    .replace("Ontology_term=", "")
                                    .split(",")
                                )
                                for hit in go:
                                    if hit not in ontologia:
                                        ontologia.append(hit)
                            if "Dbxref" in field:
                                interpro.append(
                                    field.replace('"', "").replace("Dbxref=", "")
                                )
                        if "superfamily" in db.lower():
                            superfamily_line = linha[-1].split(";")
                            for field in superfamily_line:
                                if "Name" in field:
                                    superfamily.append(
                                        field.replace('"', "").replace("Name=", "")
                                    )
                    if "\t".join(linha) == seq_reg[-1]:
                        add_none_to_empty()
                        add_to_dicts(nome_subject, ontologia, interpro, superfamily)
                        ontologia, interpro, superfamily = create_or_reset_lists()

    annotation_dict = {}

    def get_annotation(file):
        with open(file) as source:
            annotations = source.read().splitlines()
        for line in annotations:
            fields = line.split("\t", 1)
            if len(fields) != 2:
                continue
            annotation_dict[fields[0]] = re.split(r"\(InterPro|\(GO", fields[1])[
                0
            ].strip()

    transcript_dict = {}

    def get_transcripts(file):
        with open(file) as source:
            transcripts = source.read().splitlines()
        # skip the header row
        del transcripts[0]
        for line in transcripts:
            line = line.split("\t")
            transcript_dict[line[0]] = round(float(line[1]), 2)

    unique_peptide = {}
    total_peptide = {}

    def get_proteomics(file):
        with open(file) as source:
            proteomics = source.read().splitlines()
        # skip the header row
        del proteomics[0]
        for line in proteomics:
            line = line.split("\t")
            total_peptide[line[0]] = line[2]
            if line[1] != "0":
                unique_peptide[line[0]] = line[1]

    # read InterProScan evidence
    get_interpro_info(ipr_annot)
    get_interpro_info(ipr_hyp)
    get_annotation(annot)
    annotation_sources = [annotation_dict, go_dict, ipr_dict, superfamily_dict]
    # include optional analysis results
    if tr:
        get_transcripts(tr)
        annotation_sources.append(transcript_dict)
    if proteomics:
        get_proteomics(proteomics)
        annotation_sources.append(unique_peptide)
        annotation_sources.append(total_peptide)

    # dict keeps input order and gives O(1) membership checks
    unique_id_proteins = {}
    for source in annotation_sources:
        for key in source:
            unique_id_proteins[key] = None

    with open(f"AnnotaPipeline_{basename}_Summary.tsv", "w") as file_output:
        file_output.write(
            "ProteinID\tAnnotation\tIPR\tGO\tSuperfamily\tTranscript\tTPM\t"
            "Expression\tTotal_Peptide\tUnique_Peptide\n"
        )
        for protein in unique_id_proteins:
            transcript = (
                "True" if str(transcript_dict.get(protein)) != "None" else "False"
            )
            expression = (
                "True" if str(total_peptide.get(protein)) != "None" else "False"
            )
            file_output.write(
                f"{protein}\t{annotation_dict.get(protein)}\t"
                f"{str(ipr_dict.get(protein)).replace('InterPro:', '').strip()}"
                f"\t{str(go_dict.get(protein)).strip()}\t"
                f"{superfamily_dict.get(protein)}\t{transcript}"
                f"\t{transcript_dict.get(protein)}\t{expression}\t"
                f"{total_peptide.get(protein)}\t{unique_peptide.get(protein)}\n"
            )


if __name__ == "__main__":
    app()
