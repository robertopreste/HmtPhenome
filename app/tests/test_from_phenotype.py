#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import get_genes_from_phenotype, get_vars_from_phenotype, get_diseases_from_phenotype, json_from_phenotype


def test_get_genes_from_phenotype():
    expect = pd.DataFrame({
        "gene_name": ["YME1L1", "LYPLAL1", "TIMM8A", "ALDH3A2", "SLC25A4",
                      "OAT", "OTC", "ACOX1", "TRNT1", "NDUFB11", "ATAD3A",
                      "AGK", "ALDH18A1", "RECQL4", "GSTO1", "SCO2"],
        "ensembl_gene_id": ["ENSG00000136758", "ENSG00000143353",
                            "ENSG00000126953", "ENSG00000072210",
                            "ENSG00000151729", "ENSG00000065154",
                            "ENSG00000036473", "ENSG00000161533",
                            "ENSG00000072756", "ENSG00000147123",
                            "ENSG00000197785", "ENSG00000006530",
                            "ENSG00000059573", "ENSG00000160957",
                            "ENSG00000148834", "ENSG00000130489"],
        "phenotype_id": ["HP:0000545" for _ in range(16)],
        "phenotype_name": ["Myopia" for _ in range(16)]
    })
    result = get_genes_from_phenotype("HP:0000545")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_vars_from_phenotype():
    expect = pd.DataFrame({
        "chromosome": [19, 19, 19, 7, 7, 12, 12],
        "start_pos": [8373832, 8373832, 8376973, 100823256, 100823256,
                      6019472, 6019472],
        "ensembl_gene_id": ["ENSG00000167772", "ENSG00000167772",
                            "ENSG00000269386", "ENSG00000196411",
                            "ENSG00000196411", "ENSG00000110799",
                            "ENSG00000110799"],
        "gene_name": ["ANGPTL4", "ANGPTL4", "RAB11B-AS1", "EPHB4", "EPHB4",
                      "VWF", "VWF"],
        "dbsnp_id": ["rs11672433", "rs11672433", "rs1808536", "rs314308",
                     "rs314308", "rs61749397", "rs61749397"],
        "phenotype_name": ["Intracranial Hemorrhages" for _ in range(7)],
        "ref_allele": ["G", "G", "A", "C", "C", "C", "C"],
        "alt_allele": ["A", "C", "G", "G", "T", "G", "T"],
        "variant": ["chr19:8373832G>A", "chr19:8373832G>C",
                    "chr19:8376973A>G", "chr7:100823256C>G",
                    "chr7:100823256C>T", "chr12:6019472C>G",
                    "chr12:6019472C>T"],
        "phenotype_id": ["HP:0002170" for _ in range(7)]
    })
    result = get_vars_from_phenotype("HP:0002170")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_diseases_from_phenotype():
    expect = pd.DataFrame({
        "phenotype_id": ["HP:0000456" for _ in range(6)],
        "phenotype_name": ["Notched nasal tip" for _ in range(6)],
        "disease_name": ["MICROPHTHALMIA, SYNDROMIC 2; MCOPS2",
                         "ZIMMERMANN-LABAND SYNDROME 2; ZLS2",
                         "EVEN-PLUS SYNDROME; EVPLS",
                         "CRANIOFRONTONASAL SYNDROME; CFNS",
                         "FRONTONASAL DYSPLASIA 1; FND1",
                         "FRONTONASAL DYSPLASIA 2; FND2"],
        "disease_id": ["OMIM:300166", "OMIM:616455", "OMIM:616854",
                       "OMIM:304110", "OMIM:136760", "OMIM:613451"],
        "umls_disease_id": ["C1846265", "C4225321", "C4225180", "C0220767",
                            "C1876203", "C3150703"]
    })
    result = get_diseases_from_phenotype("HP:0000456")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_json_from_phenotype(): 
    expect = {"diseases": [{"disease_id": "OMIM:300166", 
                            "disease_name": "MICROPHTHALMIA, SYNDROMIC 2; MCOPS2", 
                            "phenotype_id": "HP:0000456",
                            "phenotype_name": "Notched nasal tip",
                            "umls_disease_id": "C1846265"},
                           {"disease_id": "OMIM:616455",
                            "disease_name": "ZIMMERMANN-LABAND SYNDROME 2; ZLS2",
                            "phenotype_id": "HP:0000456",
                            "phenotype_name": "Notched nasal tip",
                            "umls_disease_id": "C4225321"},
                           {"disease_id": "OMIM:616854",
                            "disease_name": "EVEN-PLUS SYNDROME; EVPLS",
                            "phenotype_id": "HP:0000456",
                            "phenotype_name": "Notched nasal tip",
                            "umls_disease_id": "C4225180"},
                           {"disease_id": "OMIM:304110",
                            "disease_name": "CRANIOFRONTONASAL SYNDROME; CFNS",
                            "phenotype_id": "HP:0000456",
                            "phenotype_name": "Notched nasal tip",
                            "umls_disease_id": "C0220767"},
                           {"disease_id": "OMIM:136760",
                            "disease_name": "FRONTONASAL DYSPLASIA 1; FND1",
                            "phenotype_id": "HP:0000456",
                            "phenotype_name": "Notched nasal tip",
                            "umls_disease_id": "C1876203"},
                           {"disease_id": "OMIM:613451",
                            "disease_name": "FRONTONASAL DYSPLASIA 2; FND2",
                            "phenotype_id": "HP:0000456",
                            "phenotype_name": "Notched nasal tip",
                            "umls_disease_id": "C3150703"}],
              "genes": [],
              "phenotype": "Notched nasal tip",
              "variants": []}
    result = json_from_phenotype("HP:0000456")
    assert result == expect
