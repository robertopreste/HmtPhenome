#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import get_gene_from_variant, get_dbsnp_from_variant, \
    get_diseases_from_dbsnp, get_phenos_from_umls, json_from_variant


def test_get_gene_from_variant():
    expect = pd.DataFrame({
        "ensembl_gene_id": ["ENSG00000198888"],
        "gene_name": ["ND1"]
    })
    result = get_gene_from_variant("M", 3308)
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_gene_from_variant_empty():
    expect = pd.DataFrame(columns=["ensembl_gene_id", "gene_name"])
    result = get_gene_from_variant("M", 17000)
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


@pytest.mark.asyncio
async def test_get_dbsnp_from_variant():
    expect = pd.DataFrame({
        "dbsnp_id": ["rs28358582", "rs28358582", "rs28358582"],
        # "variant": ["chrMT:3308T>A", "chrMT:3308T>C", "chrMT:3308T>G"]
        "variant": ["NC_012920.1:m.3308T>A", "NC_012920.1:m.3308T>C", "NC_012920.1:m.3308T>G"]
    })
    result = await get_dbsnp_from_variant("M", 3308)
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


@pytest.mark.asyncio
async def test_get_dbsnp_from_variant_empty():
    expect = pd.DataFrame(columns=["dbsnp_id", "variant"])
    result = await get_dbsnp_from_variant("M", 170000)
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_diseases_from_dbsnp():
    expect = pd.DataFrame({
        "dbsnp_id": ["rs28358582"],
        "umls_disease_id": ["C0038644"],
        "disease_id": ["OMIM:272120"],
        "disease_name": ["Sudden infant death syndrome"],
        "ass_score": [0.7]
    })
    result = get_diseases_from_dbsnp("rs28358582")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_diseases_from_dbsnp_empty():
    expect = pd.DataFrame(columns=["dbsnp_id", "umls_disease_id", "disease_id",
                                   "disease_name", "ass_score"])
    result = get_diseases_from_dbsnp("rs6666666")
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


def test_get_phenos_from_umls_empty():
    expect = pd.DataFrame(columns=["umls_disease_id", "disease_name",
                                   "disease_id", "phenotype_id",
                                   "phenotype_name"])
    result = get_phenos_from_umls("C6666666")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


@pytest.mark.asyncio
async def test_json_from_variant():
    expect = {'diseases': [{'ass_score': 0.7,
                            'dbsnp_id': 'rs28358582',
                            'disease_id': 'OMIM:272120',
                            'disease_name': 'Sudden infant death syndrome',
                            'umls_disease_id': 'C0038644'}],
              'genes': [{'ensembl_gene_id': 'ENSG00000198888',
                         'gene_name': 'ND1'}],
              'phenotypes': [],
              'variants': [{'dbsnp_id': 'rs28358582', 'disease_id': 'OMIM:272120',
                            'disease_name': 'Sudden infant death syndrome',
                            'ensembl_gene_id': 'ENSG00000198888',
                            'gene_name': 'ND1',
                            'phenotype_id': '',
                            'phenotype_name': '',
                            'umls_disease_id': 'C0038644',
                            # 'variant': 'chrMT:3308T>A'},
                            'variant': 'NC_012920.1:m.3308T>A'},
                           {'dbsnp_id': 'rs28358582', 'disease_id': 'OMIM:272120',
                            'disease_name': 'Sudden infant death syndrome',
                            'ensembl_gene_id': 'ENSG00000198888',
                            'gene_name': 'ND1',
                            'phenotype_id': '',
                            'phenotype_name': '',
                            'umls_disease_id': 'C0038644',
                            # 'variant': 'chrMT:3308T>C'},
                            'variant': 'NC_012920.1:m.3308T>C'},
                           {'dbsnp_id': 'rs28358582', 'disease_id': 'OMIM:272120',
                            'disease_name': 'Sudden infant death syndrome',
                            'ensembl_gene_id': 'ENSG00000198888',
                            'gene_name': 'ND1',
                            'phenotype_id': '',
                            'phenotype_name': '',
                            'umls_disease_id': 'C0038644',
                            # 'variant': 'chrMT:3308T>G'}]}
                            'variant': 'NC_012920.1:m.3308T>G'}]}
    result = await json_from_variant("M", 3308)
    assert result == expect
