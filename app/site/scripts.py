#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
import numpy as np
import requests
import sys
from pybiomart import Server


def get_vars_from_gene_name(gene_name):
    """Retrieve all variants associated with a specific gene, using Biomart.
    :param gene_name: name of the query gene
    :return:
    """
    # TODO: mitocarta gene list
    mitocarta_genes = pd.read_csv("")

    gene_name = gene_name.upper().lstrip("MT-")
    try:
        ens_gene_id = mitocarta_genes[mitocarta_genes["gene_symbol"] == gene_name]["ensembl_id"].values[0]
    except IndexError:  # gene not in Mitocarta
        return pd.DataFrame()

    server = Server(host="http://www.ensembl.org")
    dataset = server.marts["ENSEMBL_MART_SNP"].datasets["hsapiens_snp"]

    res = dataset.query(attributes=["chr_name", "chrom_start", "consequence_allele_string",
                                    "phenotype_description"],
                        filters={"ensembl_gene": ens_gene_id})

    res.rename({"Chromosome/scaffold name": "chromosome",
                "Chromosome/scaffold position start (bp)": "start_pos",
                "Consequence specific allele": "ref/alt allele",
                "Phenotype description": "phenotype"}, axis=1, inplace=True)
    res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
    res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
    res = res[res["allele"] != "HGMD_MUTATION"]
    res["gene"] = gene_name
    res["ensembl_gene_id"] = ens_gene_id

    res = res[["gene", "ensembl_gene_id", "chromosome", "ref_allele", "start_pos", "alt_allele",
               "phenotype"]]

    return res


def get_diseases_from_gene_name(gene_name, with_vars=False):
    """Retrieve diseases associated to a specific gene, using Ensembl.
    :param gene_name: name of the query gene
    :param with_vars: True to also retrieve phenotypes associated to variants of the gene
    :return:
    """
    gene_name = gene_name.upper().lstrip("MT-")

    server = "https://rest.ensembl.org"
    ext = "/phenotype/gene/homo_sapiens/{}?include_associated={}".format(gene_name, int(with_vars))
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        print("Wrong request.")
        sys.exit()
    res = r.json()

    df = pd.DataFrame(columns=["gene_name", "ensembl_gene_id", "disease", "phenotypes"])
    for el in res:
        row = pd.DataFrame({"gene_name": gene_name, "ensembl_gene_id": [el["Gene"]],
                            "disease": [el["description"]],
                            "phenotypes": [el["ontology_accessions"]]})
        df = df.append(row, ignore_index=True)
    df.drop_duplicates(subset="disease")

    return df


