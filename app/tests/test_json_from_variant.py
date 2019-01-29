#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import json_from_variant, get_diseases_from_variant, get_phenos_from_variant


def test_get_diseases_from_variant_m_3308():
    corr = pd.DataFrame({"ensembl_gene_id": ["ENSG00000198888", "ENSG00000198888"],
                         "gene_name": ["MT-ND1", "MT-ND1"],
                         "location": ["MT:3308-3308", "MT:3308-3308"],
                         "variation": ["rs28358582", "rs28358582"],
                         "disease": ["Carcinoma of colon", "SUDDEN INFANT DEATH SYNDROME"],
                         "phenotypes": [["EFO:0000365", "EFO:1001950"],
                                        ["EFO:0005303", "HP:0001699", "HP:0005949"]]})
    # assert get_diseases_from_variant("M", 3308).reset_index(drop=True) == corr.reset_index(drop=True)
    assert_frame_equal(get_diseases_from_variant("M", 3308).reset_index(drop=True),
                       corr.reset_index(drop=True))


def test_get_diseases_from_variant_m_420():
    corr = pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "location", "variation", "disease",
                                 "phenotypes"])
    # assert get_diseases_from_variant("M", 420) == corr
    assert_frame_equal(get_diseases_from_variant("M", 420).reset_index(drop=True),
                       corr.reset_index(drop=True))


def test_get_phenos_from_variant_m_3308():
    corr = pd.DataFrame({"chromosome": ["MT", "MT"],
                         "ref_allele": ["T", "T"],
                         "start_pos": [3308, 3308],
                         "alt_allele": ["A", "A"],
                         "phenotype": ["Carcinoma of colon", "SUDDEN INFANT DEATH SYNDROME"],
                         "phenotype_id": [["HP:0040276"], []]})
    # assert get_phenos_from_variant("M", 3308) == corr
    assert_frame_equal(get_phenos_from_variant("M", 3308).reset_index(drop=True),
                       corr.reset_index(drop=True))


def test_get_phenos_from_variant_m_420():
    corr = pd.DataFrame(columns=["chromosome", "ref_allele", "start_pos", "alt_allele", "phenotype",
                                 "phenotype_id"])
    # assert get_phenos_from_variant("M", 420) == corr
    assert_frame_equal(get_phenos_from_variant("M", 420).reset_index(drop=True),
                       corr.reset_index(drop=True))






