#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import get_vars_from_gene_id, get_diseases_from_gene  # , json_from_gene


def test_get_vars_from_gene_id():
    expect = pd.DataFrame({
        "ensembl_gene_id": ["ENSG00000210154" for _ in range(4)],
        "gene_name": ["TD" for _ in range(4)],
        "chromosome": ["MT" for _ in range(4)],
        "ref_allele": ["G", "A", "T", "T"],
        "start_pos": [7521, 7526, 7547, 7581],
        "alt_allele": ["A", "G", "C", "C"],
        "variant": ["chrMT:7521G>A", "chrMT:7526A>G",
                    "chrMT:7547T>C", "chrMT:7581T>C"],
        "dbsnp_id": ["rs200336937", "rs121434454",
                     "rs879076142", "rs201582552"]
    })
    expect.start_pos = expect.start_pos.astype("object")
    result = get_vars_from_gene_id("ENSG00000210154")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_diseases_from_gene():
    expect = pd.DataFrame({
        "ensembl_gene_id": ["ENSG00000210196"],
        "gene_name": ["TP"],
        "disease_name": ["MERRF"],
        "phenotype_ids": [["HP:0000407", "HP:0001250", "HP:0001251", "HP:0001257", "HP:0001324",
                           "HP:0001336", "HP:0001427", "HP:0002123", "HP:0002151", "HP:0003198",
                           "HP:0003200", "HP:0003542", "Orphanet:551"]]
    })
    result = get_diseases_from_gene("ENSG00000210196")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


# def test_json_from_gene():
#     # TODO
#     pass
#

