#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import json_from_gene, get_vars_from_gene_name, get_diseases_from_gene_name


def test_get_vars_from_gene_name_atp8():
    corr = pd.DataFrame({"ensembl_gene_id": ["ENSG00000228253", "ENSG00000228253",
                                             "ENSG00000228253", "ENSG00000228253"],
                         "gene_name": ["ATP8", "ATP8", "ATP8", "ATP8"],
                         "chromosome": ["MT", "MT", "MT", "MT"],
                         "ref_allele": ["C", "G", "T", "T"],
                         "start_pos": [8393, 8529, 8528, 8418],
                         "alt_allele": ["T", "A", "C", "C"],
                         "variant": ["C8393T", "G8529A", "T8528C", "T8418C"],
                         "dbsnp_id": ["rs121434446", "rs267606881", "rs387906422", "rs1057516062"],
                         "phenotype": ["Brain pseudoatrophy reversible valproate-induced susceptibility to",
                                       "Cardiomyopathy apical hypertrophic and neuropathy",
                                       "Cardiomyopathy infantile hypertrophic", "Optic neuropathy"]})
    assert_frame_equal(get_vars_from_gene_name("ATP8").reset_index(drop=True),
                       corr.reset_index(drop=True))


def test_get_diseases_from_gene_name_atp8():
    corr = pd.DataFrame({"ensembl_gene_id": ["ENSG00000228253", "ENSG00000228253", "ENSG00000228253"],
                         "gene_name": ["ATP8", "ATP8", "ATP8"],
                         "location": ["MT:8366-8572", "MT:8366-8572", "MT:8366-8572"],
                         "disease": ["Isolated ATP synthase deficiency",
                                     "Periodic paralysis with later-onset distal motor neuropathy",
                                     "KEARNS-SAYRE SYNDROME"],
                         "phenotypes": [["Orphanet:254913"], ["Orphanet:397750"], 
                                        ["HP:0000252", "HP:0000365", "HP:0000407", "HP:0000508",
                                         "HP:0000580", "HP:0000590", "HP:0000597", "HP:0000726",
                                         "HP:0000763", "HP:0000819", "HP:0000829", "HP:0000830",
                                         "HP:0001250", "HP:0001251", "HP:0001252", "HP:0001315",
                                         "HP:0001324", "HP:0001427", "HP:0001638", "HP:0001709",
                                         "HP:0001924", "HP:0001947", "HP:0001994", "HP:0002135",
                                         "HP:0002311", "HP:0002750", "HP:0002922", "HP:0003128",
                                         "HP:0003200", "HP:0003202", "HP:0003457", "HP:0004322",
                                         "HP:0004374", "HP:0007703", "HP:0008207", "HP:0011675",
                                         "Orphanet:480"]]})
    assert_frame_equal(get_diseases_from_gene_name("ATP8").reset_index(drop=True),
                       corr.reset_index(drop=True))

