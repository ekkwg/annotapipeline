#!/usr/bin/env python3
"""parse InterProScan, HMMER and RPS-BLAST evidence"""

import logging
import re
from pathlib import Path

import typer
from logging_utils import configure, redirect_stderr

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    ipr_annot: str = typer.Option(..., "-ipr_annot"),
    ipr_hyp: str = typer.Option(..., "-ipr_hyp"),
    hmm: str = typer.Option(..., "-hmm", "--hmmscan"),
    rpsblast: str = typer.Option(..., "-rpsblast"),
    basename: str = typer.Option(..., "-basename"),
) -> None:
    logger = configure(
        "FunctionalAnnotation.log",
        name="Functional Annotation",
        level=logging.DEBUG,
    )
    redirect_stderr(logger)

    def parser_interproscan(arq_entrada, arq_ipr, arq_saida):
        with open(arq_entrada) as source:
            entrada = source.read().split("##FASTA")
        interp = entrada[0].split("##sequence-region")
        del interp[0]

        db_list = ["coils", "gene3d", "mobidblite"]

        with open(arq_ipr, "a") as ipr, open(arq_saida, "a") as output:
            for seq_reg in interp:
                seq_reg = seq_reg.splitlines()
                for linha in seq_reg:
                    if linha == seq_reg[0]:
                        linha = linha.split(" ")
                        del linha[0]
                    elif linha == seq_reg[1]:
                        # skip the sequence-region header for each group
                        pass
                    else:
                        ontology = str(None)
                        anotation_db = str(None)
                        name = str(None)
                        interpro = str(None)
                        linha = linha.split("\t")
                        name_subject = linha[0]
                        db = linha[1]
                        if not any(unwanted in db.lower() for unwanted in db_list):
                            correct_db = db
                            anotation = linha[-1].split(";")
                            for query in anotation:
                                if "Ontology" in query:
                                    ontology = query.replace('"', "").replace(
                                        "Ontology_term=", ""
                                    )
                                if "signature_" in query:
                                    anotation_db = query.replace("signature_desc=", "")
                                if "Name" in query:
                                    name = query.replace('"', "").replace("Name=", "")
                                if "Dbxref" in query:
                                    interpro = query.replace('"', "").replace(
                                        "Dbxref=", ""
                                    )
                            row = (
                                f"{name_subject}\t{correct_db}\t{name}\t"
                                f"{anotation_db}\t{interpro}\t{ontology}\n"
                            )
                            ipr.write(row)
                            output.write(row)

    def pfam_format(arq_entrada, arq_saida):
        with open(arq_entrada) as source:
            lines = source.read().splitlines()
        del lines[0:3]
        del lines[-10:]
        with open(arq_saida, "w") as output:
            for line in lines:
                fields = re.sub(r"\s+", " ", line).split(" ")
                description = " ".join(fields[18:])
                del fields[18:]
                fields.append(description)
                output.write("\t".join(fields) + "\n")

    def parser_pfam(arq_entrada, arq_saida):
        with open(arq_entrada) as source:
            lines = source.read().splitlines()
        with open(arq_saida, "a") as output:
            for line in lines:
                fields = line.strip().split("\t")
                query = fields[2]
                db = "Pfam"
                ontologia = str(None)
                anotacao_db = fields[-1]
                name = fields[1]
                interpro = str(None)
                output.write(
                    f"{query}\t{db}\t{name}\t{anotacao_db}\t{interpro}\t{ontologia}\n"
                )

    def parser_rpsblast(arq_entrada, arq_rps, arq_saida):
        with open(arq_entrada) as source:
            lines = source.read().splitlines()
        with open(arq_rps, "w") as rps, open(arq_saida, "a") as output:
            for line in lines:
                hit = line.split("\t")
                query = hit[0]
                db = "CDD"
                name = hit[1]
                anotation = hit[-1]
                interpro = str(None)
                ontology = str(None)
                row = f"{query}\t{db}\t{name}\t{anotation}\t{interpro}\t{ontology}\n"
                rps.write(row)
                output.write(row)

    def sort_arq(arq_entrada, arq_saida):
        with open(arq_entrada) as source:
            lines = source.readlines()
        lines.sort()
        with open(arq_saida, "w") as ordered:
            ordered.write("ID\tDB\tDB_ACCESS\tDESCRIPTION\tIPR\tGO\n")
            ordered.writelines(str(line) for line in lines)

    # collect evidence for hypothetical proteins
    # parse InterProScan output
    parser_interproscan(
        ipr_hyp, f"InterProScan_Out_{basename}.tsv", f"Temp_{basename}.tsv"
    )
    logger.info("InterProScan parser done")
    # reformat HMMER output
    pfam_format(hmm, f"Hmmscan_Out_{basename}.tsv")
    logger.info("Hmmscan format done")
    parser_pfam(f"Hmmscan_Out_{basename}.tsv", f"Temp_{basename}.tsv")
    logger.info("Hmmscan parser done")
    # parse RPS-BLAST output
    parser_rpsblast(rpsblast, f"RPSblast_Out_{basename}.tsv", f"Temp_{basename}.tsv")
    logger.info("RPSblast parser done")
    # group and sort evidence for hypothetical proteins
    sort_arq(
        f"Temp_{basename}.tsv",
        f"{basename}_Grouped_Hypothetical_Information.tsv",
    )
    logger.info("Sorting queryes in files")

    # collect InterProScan evidence for annotated proteins
    parser_interproscan(
        ipr_annot, f"InterProScan_Out_{basename}.tsv", f"Temp_{basename}.tsv"
    )
    # remove the scratch table
    Path(f"Temp_{basename}.tsv").unlink(missing_ok=True)
    logger.info("Temporary file removed")
    logger.info("Functional Annotation step is completed")


if __name__ == "__main__":
    app()
