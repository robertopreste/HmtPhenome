#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
from app.site.scripts import create_variant_string, pheno_name_to_id, \
    pheno_id_to_term, ensembl_gene_id_to_entrez, disease_id_to_name, \
    disease_name_to_id, create_variant_HGVS


def test_mt_snp_variant():
    chrom = "M"
    nt_start = 3308
    ref_all = "G"
    alt_all = "A"
    expect = "chr{}:{}{}>{}".format(chrom, nt_start, ref_all, alt_all)
    result = create_variant_string(chrom, nt_start, ref_all, alt_all)

    assert result == expect


def test_mt_del_variant():
    chrom = "M"
    nt_start = 3308
    ref_all = "G"
    alt_all = "-"
    expect = "chr{}:{}del{}".format(chrom, nt_start, ref_all)
    result = create_variant_string(chrom, nt_start, ref_all, alt_all)

    assert result == expect


def test_mt_ins_variant():
    chrom = "M"
    nt_start = 3308
    ref_all = "-"
    alt_all = "A"
    expect = "chr{}:{}ins{}".format(chrom, nt_start, alt_all)
    result = create_variant_string(chrom, nt_start, ref_all, alt_all)

    assert result == expect


def test_mt_snp_variant_HGVS():
    chrom = "M"
    nt_start = 3308
    ref_all = "G"
    alt_all = "A"
    expect = "NC_012920.1:m.3308G>A"
    result = create_variant_HGVS(chrom, nt_start, ref_all, alt_all)

    assert result == expect


def test_mt_del_variant_HGVS():
    chrom = "M"
    nt_start = 3308
    ref_all = "G"
    alt_all = "-"
    expect = "NC_012920.1:m.3308_3309del"
    result = create_variant_HGVS(chrom, nt_start, ref_all, alt_all)

    assert result == expect


def test_mt_ins_variant_HGVS():
    chrom = "M"
    nt_start = 3308
    ref_all = "-"
    alt_all = "A"
    expect = "NC_012920.1:m.3308_3309insA"
    result = create_variant_HGVS(chrom, nt_start, ref_all, alt_all)

    assert result == expect


@pytest.mark.asyncio
async def test_pheno_name_to_id():
    pheno_name = "Nausea"
    expect = ["HP:0002017", "HP:0002018"]
    result = await pheno_name_to_id(pheno_name)

    assert result == expect


@pytest.mark.asyncio
async def test_pheno_id_to_term_HP():
    pheno_id = "HP:0000545"
    expect = "Myopia"
    result = await pheno_id_to_term(pheno_id)

    assert result == expect


def test_pheno_id_to_term_Orphanet():
    pass


def test_pheno_id_to_term_EFO():
    pass


def test_disease_id_to_name_OMIM():
    expect = "ADULT SYNDROME"
    result = disease_id_to_name("OMIM:103285")

    assert result == expect


def test_disease_id_to_name_Orpha():
    expect = "Oligodontia"
    result = disease_id_to_name("ORPHA:99798")

    assert result == expect


def test_disease_name_to_id_OMIM():
    expect = "OMIM:103285"
    result = disease_name_to_id("ADULT SYNDROME")

    assert result == expect


def test_disease_name_to_id_Orpha():
    expect = "ORPHA:99798"
    result = disease_name_to_id("Oligodontia")

    assert result == expect


@pytest.mark.asyncio
async def test_ensembl_gene_id_to_entrez():
    expect = pd.DataFrame({
        "ensembl_gene_id": ["ENSG00000198888"],
        "gene_name": ["MT-ND1"],
        "entrez_gene_id": [4535]
    })
    result = await ensembl_gene_id_to_entrez("ENSG00000198888")

    assert_frame_equal(result.reset_index(drop=True),
                       expect.reset_index(drop=True))

