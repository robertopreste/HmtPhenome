#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import get_vars_from_gene, get_diseases_from_gene, \
    json_from_gene


def test_get_vars_from_gene():
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
    result = get_vars_from_gene("ENSG00000210154")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_vars_from_gene_empty():
    expect = pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                   "chromosome", "ref_allele", "start_pos",
                                   "alt_allele", "variant", "dbsnp_id"])
    result = get_vars_from_gene("ENSG6666666")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_diseases_from_gene():
    expect = pd.DataFrame({
        "ensembl_gene_id": ["ENSG00000210196"],
        "gene_name": ["TP"],
        "disease_name": ["MERRF"],
        "phenotype_ids": [["HP:0000407", "HP:0001250", "HP:0001251",
                           "HP:0001257", "HP:0001324", "HP:0001336",
                           "HP:0001427", "HP:0002123", "HP:0002151",
                           "HP:0003198", "HP:0003200", "HP:0003542",
                           "Orphanet:551"]]
    })
    result = get_diseases_from_gene("ENSG00000210196")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_get_diseases_from_gene_empty():
    expect = pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                   "disease_name", "phenotype_ids"])
    result = get_diseases_from_gene("ENSG66666666")
    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))


def test_json_from_gene():
    expect = {'diseases': [{'disease_id': 'ORPHA:551', 'disease_name': 'MERRF',
                            'ensembl_gene_id': 'ENSG00000210196',
                            'gene_name': 'TP',
                            'phenotype_ids': ['HP:0000407', 'HP:0001250',
                                              'HP:0001251', 'HP:0001257',
                                              'HP:0001324', 'HP:0001336',
                                              'HP:0001427', 'HP:0002123',
                                              'HP:0002151', 'HP:0003198',
                                              'HP:0003200', 'HP:0003542',
                                              'Orphanet:551'],
                            'phenotype_names': ['Sensorineural hearing impairment',
                                                'Seizures', 'Ataxia', 'Spasticity',
                                                'Muscle weakness', 'Myoclonus',
                                                'Mitochondrial inheritance',
                                                'Generalized myoclonic seizures',
                                                'Increased serum lactate',
                                                'Myopathy', 'Ragged-red muscle fibers',
                                                'Increased serum pyruvate', 'MERRF'],
                            'umls_disease_id': ''}],
              'genes': [{'ensembl_gene_id': 'ENSG00000210196',
                         'gene_name': 'TP'}],
              'phenotypes': [{'phenotype_id': 'HP:0000407',
                              'phenotype_name': 'Sensorineural hearing impairment'},
                             {'phenotype_id': 'HP:0001250',
                              'phenotype_name': 'Seizures'},
                             {'phenotype_id': 'HP:0001251',
                              'phenotype_name': 'Ataxia'},
                             {'phenotype_id': 'HP:0001257',
                              'phenotype_name': 'Spasticity'},
                             {'phenotype_id': 'HP:0001324',
                              'phenotype_name': 'Muscle weakness'},
                             {'phenotype_id': 'HP:0001336',
                              'phenotype_name': 'Myoclonus'},
                             {'phenotype_id': 'HP:0001427',
                              'phenotype_name': 'Mitochondrial inheritance'},
                             {'phenotype_id': 'HP:0002123',
                              'phenotype_name': 'Generalized myoclonic seizures'},
                             {'phenotype_id': 'HP:0002151',
                              'phenotype_name': 'Increased serum lactate'},
                             {'phenotype_id': 'HP:0003198',
                              'phenotype_name': 'Myopathy'},
                             {'phenotype_id': 'HP:0003200',
                              'phenotype_name': 'Ragged-red muscle fibers'},
                             {'phenotype_id': 'HP:0003542',
                              'phenotype_name': 'Increased serum pyruvate'},
                             {'phenotype_id': 'Orphanet:551',
                              'phenotype_name': 'MERRF'}],
              'variants': [{'chromosome': 'MT', 'dbsnp_id': 'rs199474700',
                            'disease_name': 'PARKINSON DISEASE, MITOCHONDRIAL (disorder)',
                            'ensembl_gene_id': 'ENSG00000210196',
                            'gene_name': 'TP',
                            'umls_disease_id': 'C1838867',
                            'variant': 'chrMT:15965A>G'},
                           {'chromosome': 'MT', 'dbsnp_id': 'rs199474701',
                            'disease_name': 'MERFF SYNDROME',
                            'ensembl_gene_id': 'ENSG00000210196',
                            'gene_name': 'TP', 'umls_disease_id': 'C4016625',
                            'variant': 'chrMT:15967G>A'},
                           {'chromosome': 'MT', 'dbsnp_id': 'rs375213730',
                            'disease_name': '',
                            'ensembl_gene_id': 'ENSG00000210196',
                            'gene_name': 'TP', 'umls_disease_id': '',
                            'variant': 'chrMT:15970T>C'},
                           {'chromosome': 'MT', 'dbsnp_id':
                               'rs201041059', 'disease_name': '',
                            'ensembl_gene_id': 'ENSG00000210196',
                            'gene_name': 'TP',
                            'umls_disease_id': '', 'variant': 'chrMT:15978C>T'},
                           {'chromosome': 'MT', 'dbsnp_id': 'rs199474699',
                            'disease_name': 'Myopathy',
                            'ensembl_gene_id': 'ENSG00000210196',
                            'gene_name': 'TP', 'umls_disease_id': 'C0026848',
                            'variant': 'chrMT:15990C>T'},
                           {'chromosome': 'MT', 'dbsnp_id': 'rs201864830',
                            'disease_name': '',
                            'ensembl_gene_id': 'ENSG00000210196',
                            'gene_name': 'TP', 'umls_disease_id': '',
                            'variant': 'chrMT:16017T>C'},
                           {'chromosome': 'MT', 'dbsnp_id': 'rs55934780',
                            'disease_name': '',
                            'ensembl_gene_id': 'ENSG00000210196',
                            'gene_name': 'TP', 'umls_disease_id': '',
                            'variant': 'chrMT:16023G>A'}]}
    result = json_from_gene("ENSG00000210196")
    assert result == expect
