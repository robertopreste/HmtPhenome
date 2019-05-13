#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import get_umls_from_disease_id, \
    get_genes_from_disease_id, get_vars_from_disease_id, \
    get_phenos_from_disease_id, json_from_disease


def test_get_umls_from_disease_id():
    expect = "C3150703"
    result = get_umls_from_disease_id("OMIM:613451")
    assert result == expect


# def test_get_genes_from_disease_id_OMIM():
#     expect = pd.DataFrame({
#         "ass_score": [0.0, 0.0],
#         "disease_id": ["OMIM:101400", "OMIM:101400"],
#         "disease_name": ["SAETHRE-CHOTZEN SYNDROME; SCS",
#                          "SAETHRE-CHOTZEN SYNDROME; SCS"],
#         "ensembl_gene_id": ["", ""],
#         "entrez_gene_id": ["2263", "7291"],
#         "gene_name": ["FGFR2", "TWIST1"],
#         "umls_disease_id": ["C1863371", "C1863371"]
#     })
#     result = get_genes_from_disease_id("OMIM:101400")
#     assert_frame_equal(result.reset_index(drop=True),
#                        expect.reset_index(drop=True))
#
#
# def test_get_genes_from_disease_id_Orpha():
#     expect = pd.DataFrame({
#         "ass_score": [0.0],
#         "disease_id": ["ORPHA:239"],
#         "disease_name": ["Dyggve-Melchior-Clausen disease"],
#         "ensembl_gene_id": [""],
#         "entrez_gene_id": ["54808"],
#         "gene_name": ["DYM"],
#         "umls_disease_id": [""]
#     })
#     result = get_genes_from_disease_id("ORPHA:239")
#     assert_frame_equal(result.reset_index(drop=True),
#                        expect.reset_index(drop=True))


def test_get_genes_from_disease_id_empty():
    expect = pd.DataFrame(columns=["umls_disease_id", "disease_name",
                                   "disease_id", "entrez_gene_id", "gene_name",
                                   "ensembl_gene_id", "ass_score"])
    result = get_genes_from_disease_id("OMIM:666666666")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_vars_from_disease_id():
    expect = pd.DataFrame({
        "ensembl_gene_id": ["ENSG00000122691" for _ in range(3)],
        "gene_name": ["TWIST1" for _ in range(3)],
        "dbsnp_id": ["rs104894057", "rs104894059", "rs121909189"],
        "variant": ["chr7:19116966T>G", "chr7:19116856T>C", "chr7:19116930A>G"],
        "umls_disease_id": ["C1863371" for _ in range(3)],
        "disease_name": ["BLEPHAROPHIMOSIS, EPICANTHUS INVERSUS, AND PTOSIS 3, FORMERLY"
                         for _ in range(3)],
        "disease_id": ["OMIM:101400" for _ in range(3)]
    })
    result = get_vars_from_disease_id("OMIM:101400")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_vars_from_disease_id_empty():
    expect = pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "dbsnp_id",
                                   "variant", "umls_disease_id",
                                   "disease_name", "disease_id"])
    result = get_vars_from_disease_id("OMIM:66666666")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_phenos_from_disease_id():
    expect = pd.DataFrame({
        "umls_disease_id": ["C3150703" for _ in range(22)],
        "disease_name": ["FRONTONASAL DYSPLASIA 2" for _ in range(22)],
        "disease_id": ["OMIM:613451" for _ in range(22)],
        "phenotype_id": ["HP:0000457", "HP:0003828", "HP:0001320",
                         "HP:0000633", "HP:0000582", "HP:0002697",
                         "HP:0010761", "HP:0012745", "HP:0000316",
                         "HP:0000653", "HP:0001249", "HP:0001596",
                         "HP:0000007", "HP:0002084", "HP:0000456",
                         "HP:0000252", "HP:0002079", "HP:0000431",
                         "HP:0000535", "HP:0005280", "HP:0001363",
                         "HP:0000966"],
        "phenotype_name": ["Depressed nasal ridge", "Variable expressivity",
                           "Cerebellar vermis hypoplasia",
                           "Decreased lacrimation",
                           "Upslanted palpebral fissure", "Parietal foramina",
                           "Broad columella", "Short palpebral fissure",
                           "Hypertelorism", "Sparse eyelashes",
                           "Intellectual disability", "Alopecia",
                           "Autosomal recessive inheritance", "Encephalocele",
                           "Bifid nasal tip", "Microcephaly",
                           "Hypoplasia of the corpus callosum",
                           "Wide nasal bridge", "Sparse and thin eyebrow",
                           "Depressed nasal bridge", "Craniosynostosis",
                           "Hypohidrosis"]
    })
    result = get_phenos_from_disease_id("OMIM:613451")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_phenos_from_disease_id_empty():
    expect = pd.DataFrame(columns=["umls_disease_id", "disease_name",
                                   "disease_id", "phenotype_id",
                                   "phenotype_name"])
    result = get_phenos_from_disease_id("OMIM:666666666")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))
