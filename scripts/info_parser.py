#!/usr/bin/env python3
"""merge annotation sources into the product table"""

import os
import os.path
from pathlib import Path

import typer

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "-help", "--help"]},
)


@app.command()
def main(
    ipr1: str = typer.Option(..., "-ipr1"),
    ipr2: str = typer.Option(..., "-ipr2"),
    annot: str = typer.Option(..., "-a", "--annotated"),
    nohit: str = typer.Option(..., "-nh", "--no_hit"),
    hypo: str = typer.Option(..., "-hy", "--hypothetical"),
) -> None:
    hypo_ids: list[str] = []
    nohit_ids: list[str] = []

    def temporary_query(arq):
        # append a sentinel BLAST row so the last query is processed
        try:
            temp = arq[-1].split("\t")
            del temp[0]
            temp.insert(0, "QueryTemp")
            arq.append("\t".join(temp))
        # an empty InterProScan result has no row to copy
        except IndexError:
            arq.append("\t".join("QueryTemp"))

    def parser_interproscan(arq_entrada, arq_ipr):
        with open(arq_entrada) as source:
            entrada = source.read().split("##FASTA")
        interp = entrada[0].split("##sequence-region")
        del interp[0]

        lista = ["coils", "gene3d", "mobidblite"]

        with open(arq_ipr, "a") as ipr:
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
                        ontologia = str(None)
                        anotacao_db = str(None)
                        name = str(None)
                        interpro = str(None)
                        linha = linha.split("\t")
                        nome_subject = linha[0]
                        db = linha[1]
                        try:
                            evalue = str(linha[5]).replace(",", ".").lower()
                        except Exception:
                            # keep the source format when normalisation fails
                            evalue = str(linha[5])
                        if not any(s in db.lower() for s in lista):
                            db_certo = db
                            anotacao = linha[-1].split(";")
                            for a in anotacao:
                                if "Ontology" in a:
                                    ontologia = a.replace('"', "").replace(
                                        "Ontology_term=", ""
                                    )
                                if "signature_" in a:
                                    anotacao_db = a.replace("signature_desc=", "")
                                if "Name" in a:
                                    name = a.replace('"', "").replace("Name=", "")
                                if "Dbxref" in a:
                                    interpro = a.replace('"', "").replace("Dbxref=", "")
                            ipr.write(
                                f"{nome_subject}\t{db_certo}\t{name}\t"
                                f"{evalue}\t{anotacao_db}\t{interpro}\t"
                                f"{ontologia}\n"
                            )

    # write products without InterPro or GO evidence
    def write_no_ipr(ids_dict):
        # annotated products without an InterProScan hit keep their description
        if len(ids_dict.keys()) > 0:
            for anot in ids_dict:
                output.write(f"{anot}\t{str(ids_dict.get(anot)).strip()}\n")
        if len(hypo_ids) > 0:
            output.writelines(f"{hyp}\thypothetical protein\n" for hyp in hypo_ids)
        if len(nohit_ids) > 0:
            output.writelines(f"{hyp}\thypothetical protein\n" for hyp in nohit_ids)

    # attach InterPro and GO terms to hypothetical or no-hit products
    def intepro_process(ids_dict):
        # read the preprocessed identifier, IPR and GO table
        with open("Interpro_out_tmp.txt") as source:
            interpro_out = source.read().splitlines()
        temporary_query(interpro_out)

        old_id = interpro_out[0].split("\t")
        old_id = old_id[0]

        # collect terms for the current identifier
        go_list = []
        ipr_list = []
        for query in interpro_out:
            title = query.split("\t")
            new_id = title[0]
            ipr = str(title[5]).split(",")
            go = str(title[6]).split(",")
            if old_id != new_id:
                # sort terms by their numeric suffix
                try:
                    ipr_list.sort(key=lambda item: str(item).split("IPR")[1])
                    go_list.sort(key=lambda item: str(item).split("GO")[1])
                except Exception:
                    # missing terms are expected, not a warning
                    pass
                # annotated products keep their BLAST description
                if old_id in ids_dict:
                    if (len(ipr_list) > 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({','.join(ipr_list)!s},{','.join(go_list)!s})\n"
                        )
                    elif (len(ipr_list) > 0) and (len(go_list) == 0):
                        output.write(
                            f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({','.join(ipr_list)!s})\n"
                        )
                    elif (len(ipr_list) == 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\t{str(ids_dict.get(old_id)).strip()} ({','.join(go_list)!s})\n"
                        )
                    else:
                        output.write(f"{old_id}\t{str(ids_dict.get(old_id)).strip()}\n")
                    del ids_dict[old_id]
                # otherwise treat the product as hypothetical
                elif old_id in hypo_ids:
                    if (len(ipr_list) > 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(ipr_list)!s},{','.join(go_list)!s})\n"
                        )
                    elif (len(ipr_list) > 0) and (len(go_list) == 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(ipr_list)!s})\n"
                        )
                    elif (len(ipr_list) == 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(go_list)!s})\n"
                        )
                    else:
                        output.write(f"{old_id}\thypothetical protein\n")
                    hypo_ids.remove(old_id)
                else:
                    if (len(ipr_list) > 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(ipr_list)!s},{','.join(go_list)!s})\n"
                        )
                    elif (len(ipr_list) > 0) and (len(go_list) == 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(ipr_list)!s})\n"
                        )
                    elif (len(ipr_list) == 0) and (len(go_list) > 0):
                        output.write(
                            f"{old_id}\thypothetical protein ({','.join(go_list)!s})\n"
                        )
                    else:
                        output.write(f"{old_id}\thypothetical protein\n")
                    nohit_ids.remove(old_id)
                go_list.clear()
                ipr_list.clear()
                # seed the term lists for the next identifier
                for match in ipr:
                    if (match != "None") and (match not in ipr_list):
                        ipr_list.append(match)
                for match in go:
                    if (match != "None") and (match not in go_list):
                        go_list.append(match)
            else:
                for match in ipr:
                    if (match != "None") and (match not in ipr_list):
                        ipr_list.append(match)
                for match in go:
                    if (match != "None") and (match not in go_list):
                        go_list.append(match)
            old_id = new_id
        write_no_ipr(ids_dict)
        os.remove("Interpro_out_tmp.txt")

    Path("Interpro_out_tmp.txt").unlink(missing_ok=True)
    # build the preprocessed InterProScan table
    if (os.path.getsize(ipr1) == 0) and (os.path.getsize(ipr2) == 0):
        # skip parsing when both files are empty
        pass
    else:
        # parse only files with results
        if os.path.getsize(ipr1) != 0:
            parser_interproscan(ipr1, "Interpro_out_tmp.txt")
        if os.path.getsize(ipr2) != 0:
            parser_interproscan(ipr2, "Interpro_out_tmp.txt")

    # read annotated products into an identifier map
    with open(annot) as source:
        annotation_lines = source.read().splitlines()
    ids = []
    desc = []
    for anot in annotation_lines:
        anot = anot.split("\t", 1)
        if len(anot) != 2:
            continue
        ids.append(anot[0])
        desc.append(anot[1])

    ids_dict = dict(
        zip(ids, desc, strict=False)
    )  # map each identifier to its description

    # read hypothetical identifiers
    if os.path.getsize(hypo) != 0:
        with open(hypo) as source:
            hypo_ids = source.read().splitlines()
    # read no-hit identifiers
    if os.path.getsize(nohit) != 0:
        with open(nohit) as source:
            nohit_ids = source.read().splitlines()

    with open("All_annotation_products.txt", "w") as output:
        if os.path.isfile("Interpro_out_tmp.txt") and os.path.getsize(
            "Interpro_out_tmp.txt"
        ):
            intepro_process(ids_dict)
        else:
            write_no_ipr(ids_dict)
            if os.path.isfile("Interpro_out_tmp.txt"):
                os.remove("Interpro_out_tmp.txt")


if __name__ == "__main__":
    app()
