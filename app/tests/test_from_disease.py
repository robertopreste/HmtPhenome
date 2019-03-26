#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import get_umls_from_disease_id, get_genes_from_disease_id


def test_get_umls_from_disease_id():
    expect = "C3150703"
    result = get_umls_from_disease_id("OMIM:613451")
    assert result == expect


def test_get_genes_from_disease_id_OMIM():
    expect = pd.DataFrame({
        "ass_score": [0.0, 0.0],
        "disease_id": ["OMIM:101400", "OMIM:101400"],
        "disease_name": ["SAETHRE-CHOTZEN SYNDROME; SCS",
                         "SAETHRE-CHOTZEN SYNDROME; SCS"],
        "ensembl_gene_id": ["", ""],
        "entrez_gene_id": ["2263", "7291"],
        "gene_name": ["FGFR2", "TWIST1"],
        "umls_disease_id": ["C1863371", "C1863371"]
    })
    result = get_genes_from_disease_id("OMIM:101400")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_genes_from_disease_id_Orpha():
    expect = pd.DataFrame({
        "ass_score": [0.0],
        "disease_id": ["ORPHA:239"],
        "disease_name": ["Dyggve-Melchior-Clausen disease"],
        "ensembl_gene_id": [""],
        "entrez_gene_id": ["54808"],
        "gene_name": ["DYM"],
        "umls_disease_id": [""]
    })
    result = get_genes_from_disease_id("ORPHA:239")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))







