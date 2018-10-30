#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from .models import Mitocarta
# import pandas as pd
# import numpy as np
# import requests
# import sys
# from pybiomart import Server


def get_genes():
    chr_dict = {}
    chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
            "chr20", "chr21", "chr22", "chrX", "chrM"]
    for chrom in chrs:
        lista = set()
        q = Mitocarta.query.filter(Mitocarta.hg_chr == chrom).all()
        for el in q:
            lista.add(el.gene_symbol)
        chr_dict[chrom] = sorted(list(lista))

    return chr_dict


# def get_vars_from_gene_name(gene_name):
#     """Retrieve all variants associated with a specific gene, using Biomart.
#     :param gene_name: name of the query gene
#     :return:
#     """
#     # TODO: mitocarta gene list
#     mitocarta_genes = pd.read_csv("")
#
#     gene_name = gene_name.upper().lstrip("MT-")
#     try:
#         ens_gene_id = mitocarta_genes[mitocarta_genes["gene_symbol"] == gene_name]["ensembl_id"].values[0]
#     except IndexError:  # gene not in Mitocarta
#         return pd.DataFrame()
#
#     server = Server(host="http://www.ensembl.org")
#     dataset = server.marts["ENSEMBL_MART_SNP"].datasets["hsapiens_snp"]
#
#     res = dataset.query(attributes=["chr_name", "chrom_start", "consequence_allele_string",
#                                     "phenotype_description"],
#                         filters={"ensembl_gene": ens_gene_id})
#
#     res.rename({"Chromosome/scaffold name": "chromosome",
#                 "Chromosome/scaffold position start (bp)": "start_pos",
#                 "Consequence specific allele": "ref/alt allele",
#                 "Phenotype description": "phenotype"}, axis=1, inplace=True)
#     res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
#     res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
#     res = res[res["allele"] != "HGMD_MUTATION"]
#     res["gene"] = gene_name
#     res["ensembl_gene_id"] = ens_gene_id
#
#     res = res[["gene", "ensembl_gene_id", "chromosome", "ref_allele", "start_pos", "alt_allele",
#                "phenotype"]]
#
#     return res
#
#
# def get_diseases_from_gene_name(gene_name, with_vars=False):
#     """Retrieve diseases associated to a specific gene, using Ensembl.
#     :param gene_name: name of the query gene
#     :param with_vars: True to also retrieve phenotypes associated to variants of the gene
#     :return:
#     """
#     gene_name = gene_name.upper().lstrip("MT-")
#
#     server = "https://rest.ensembl.org"
#     ext = "/phenotype/gene/homo_sapiens/{}?include_associated={}".format(gene_name, int(with_vars))
#     r = requests.get(server + ext, headers={"Content-Type": "application/json"})
#
#     if not r.ok:
#         r.raise_for_status()
#         print("Wrong request.")
#         sys.exit()
#     res = r.json()
#
#     df = pd.DataFrame(columns=["gene_name", "ensembl_gene_id", "location", "disease", "phenotypes"])
#     for el in res:
#         row = pd.DataFrame({"gene_name": gene_name, "ensembl_gene_id": [el["Gene"]],
#                             "location": el["location"], "disease": [el["description"]],
#                             "phenotypes": [el["ontology_accessions"]]})
#         df = df.append(row, ignore_index=True)
#     df.drop_duplicates(subset="disease", inplace=True)
#
#     return df
#
#
# def get_gene_from_variant(chrom, var_start, var_end=None):
#     """Retrieve the gene to which the provided variant belongs, using Biomart.
#     :param chrom: chromosome name (1:22, X, Y, MT)
#     :param var_start: variant starting position
#     :param var_end: variant ending position
#     :return:
#     """
#
#     chrom = chrom.upper()
#     if var_end is None:
#         var_end = var_start
#
#     server = Server(host="http://www.ensembl.org")
#     dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]
#
#     res = dataset.query(attributes=["ensembl_gene_id", "external_gene_name"],
#                         filters={"chromosome_name": chrom, "start": var_start, "end": var_end})
#
#     res.rename({"Gene stable ID": "ensembl_gene_id", "Gene name": "gene_name"},
#                axis=1, inplace=True)
#
#     return res
#
#
# def get_pheno_from_variant(chrom, var_start, var_end=None):
#     """Retrieve phenotypes associated with a specific variant, using Biomart.
#     :param chrom: chromosome name (1:22, X, Y, MT)
#     :param var_start: variant starting position
#     :param var_end: variant ending position
#     :return:
#     """
#
#     chrom = chrom.upper()
#     if var_end is None:
#         var_end = var_start
#
#     server = Server(host="http://www.ensembl.org")
#     dataset = server.marts["ENSEMBL_MART_SNP"].datasets["hsapiens_snp"]
#
#     res = dataset.query(attributes=["chr_name", "chrom_start", "consequence_allele_string",
#                                     "phenotype_description"],
#                         filters={"chr_name": chrom, "start": var_start, "end": var_end})
#
#     res.rename({"Chromosome/scaffold name": "chromosome",
#                 "Chromosome/scaffold position start (bp)": "start_pos",
#                 "Consequence specific allele": "ref/alt allele",
#                 "Phenotype description": "phenotype"}, axis=1, inplace=True)
#     res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
#     res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
#     res = res[res["alt_allele"] != "HGMD_MUTATION"]
#     res = res[["chromosome", "ref_allele", "start_pos", "alt_allele", "phenotype"]]
#     res.drop_duplicates(subset="phenotype", inplace=True)
#
#     return res
#
#
# def get_genes_from_disease_name(disease_name):
#     """Retrieve genes related to a disease, using Biomart.
#     :param disease_name: name of the query disease
#     :return:
#     """
#
#     server = Server(host="http://www.ensembl.org")
#     dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]
#     res = dataset.query(attributes=["ensembl_gene_id", "external_gene_name", "description",
#                                     "chromosome_name", "start_position", "end_position",
#                                     "phenotype_description"],
#                         filters={"phenotype_description": disease_name})
#
#     res.rename({"Gene stable ID": "ensembl_gene_id", "Gene name": "gene_name",
#                 "Gene description": "gene_descr", "Chromosome/scaffold name": "chromosome",
#                 "Gene start (bp)": "start_pos", "Gene end (bp)": "stop_pos",
#                 "Phenotype description": "phenotype"}, axis=1, inplace=True)
#     res.drop_duplicates(inplace=True)
#
#     return res




