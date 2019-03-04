#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import get_gene_from_variant, get_dbsnp_from_variant, get_diseases_from_dbsnp, get_phenos_from_umls, json_from_variant, network_from_variant_json


def test_get_gene_from_variant():
    expect = pd.DataFrame({
        "ensembl_gene_id": ["ENSG00000198888"],
        "gene_name": ["ND1"]
    })
    result = get_gene_from_variant("M", 3308)
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_dbsnp_from_variant():
    expect = pd.DataFrame({
        "dbsnp_id": ["rs28358582", "rs28358582", "rs28358582"],
        "variant": ["chrMT:3308T>A", "chrMT:3308T>C", "chrMT:3308T>G"]
    })
    result = get_dbsnp_from_variant("M", 3308)
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_diseases_from_dbsnp():
    expect = pd.DataFrame({
        "dbsnp_id": ["rs28358582"],
        "umls_disease_id": ["C0038644"],
        "disease_name": ["Sudden infant death syndrome"],
        "ass_score": [0.7]
    })
    result = get_diseases_from_dbsnp("rs28358582")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_phenos_from_umls():
    expect = pd.DataFrame({
        "umls_disease_id": ["C3665349" for _ in range(14)],
        "disease_name": ["Secondary hypothyroidism" for _ in range(14)],
        "disease_id": ["OMIM:275100" for _ in range(14)],
        "phenotype_id": ["HP:0000007", "HP:0001252", "HP:0005280", "HP:0001537", "HP:0000158",
                         "HP:0001290", "HP:0001615", "HP:0010864", "HP:0000851", "HP:0006887",
                         "HP:0008850", "HP:0000260", "HP:0001539", "HP:0001939"],
        "phenotype_name": ["Autosomal recessive inheritance", "Muscular hypotonia",
                           "Depressed nasal bridge", "Umbilical hernia", "Macroglossia",
                           "Generalized hypotonia", "Hoarse cry", "Intellectual disability, severe",
                           "Congenital hypothyroidism", "Intellectual disability, progressive",
                           "Severe postnatal growth retardation", "Wide anterior fontanel",
                           "Omphalocele", "Abnormality of metabolism/homeostasis"]
    })
    result = get_phenos_from_umls("C3665349")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_json_from_variant():
    expect = {'diseases': [{'ass_score': 0.7,
                            'dbsnp_id': 'rs28358582',
                            'disease_name': 'Sudden infant death syndrome',
                            'umls_disease_id': 'C0038644'}],
              'genes': [{'ensembl_gene_id': 'ENSG00000198888',
                         'gene_name': 'ND1'}],
              'phenotypes': [],
              'variants': [{'dbsnp_id': 'rs28358582',
                            'disease_name': 'Sudden infant death syndrome',
                            'ensembl_gene_id': 'ENSG00000198888',
                            'gene_name': 'ND1',
                            'phenotype_id': '',
                            'phenotype_name': '',
                            'umls_disease_id': 'C0038644',
                            'variant': 'chrMT:3308T>A'},
                           {'dbsnp_id': 'rs28358582',
                            'disease_name': 'Sudden infant death syndrome',
                            'ensembl_gene_id': 'ENSG00000198888',
                            'gene_name': 'ND1',
                            'phenotype_id': '',
                            'phenotype_name': '',
                            'umls_disease_id': 'C0038644',
                            'variant': 'chrMT:3308T>C'},
                           {'dbsnp_id': 'rs28358582',
                            'disease_name': 'Sudden infant death syndrome',
                            'ensembl_gene_id': 'ENSG00000198888',
                            'gene_name': 'ND1',
                            'phenotype_id': '',
                            'phenotype_name': '',
                            'umls_disease_id': 'C0038644',
                            'variant': 'chrMT:3308T>G'}]}
    result = json_from_variant("M", 3308)
    assert result == expect


# TODO: some items are sorted differently in the results
# def test_network_from_variant_json():
#     expect = {'diseases': {('Sudden infant death syndrome', 'C0038644')},
#               'edges': [{'from': 3, 'to': 4},
#                         {'from': 3, 'to': 5},
#                         {'from': 2, 'to': 4},
#                         {'from': 2, 'to': 5},
#                         {'from': 1, 'to': 4},
#                         {'from': 1, 'to': 5}],
#               'genes': {('ND1', 'ENSG00000198888')},
#               'nodes': [{'color': {'background': '#F9CF45', 'border': '#CCAA39'},
#                          'id': 1, 'label': 'chrMT:3308T>G'},
#                         {'color': {'background': '#F9CF45', 'border': '#CCAA39'},
#                          'id': 2, 'label': 'chrMT:3308T>C'},
#                         {'color': {'background': '#F9CF45', 'border': '#CCAA39'},
#                          'id': 3, 'label': 'chrMT:3308T>A'},
#                         {'color': {'background': '#739E82', 'border': '#5F826B'},
#                          'id': 4, 'label': 'ND1'},
#                         {'color': {'background': '#D7816A', 'border': '#B06A57'},
#                          'id': 5, 'label': 'Sudden infant death syndrome'}],
#               'phenotypes': set(),
#               'variants': {('chrMT:3308T>A', 'rs28358582', 'ND1'),
#                            ('chrMT:3308T>C', 'rs28358582', 'ND1'),
#                            ('chrMT:3308T>G', 'rs28358582', 'ND1')}}
#     result = network_from_variant_json(json_from_variant("M", 3308))
#     assert result == expect

