#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pytest
from app.site.scripts import create_variant_string, pheno_name_to_id, pheno_id_to_term


def test_mt_snp_variant():
    chrom = "M"
    nt_start = 3308
    ref_all = "G"
    alt_all = "A"
    expect = "chr{}:{}{}>{}".format(chrom, nt_start, ref_all, alt_all)

    assert create_variant_string(chrom, nt_start, ref_all, alt_all) == expect


def test_mt_del_variant():
    chrom = "M"
    nt_start = 3308
    ref_all = "G"
    alt_all = "-"
    expect = "chr{}:{}del{}".format(chrom, nt_start, ref_all)

    assert create_variant_string(chrom, nt_start, ref_all, alt_all) == expect


def test_mt_ins_variant():
    chrom = "M"
    nt_start = 3308
    ref_all = "-"
    alt_all = "A"
    expect = "chr{}:{}ins{}".format(chrom, nt_start, alt_all)

    assert create_variant_string(chrom, nt_start, ref_all, alt_all) == expect


def test_pheno_name_to_id():
    pheno_name = "Myopia"
    expect = ["HP:0000545", "HP:0031624", "HP:0025573", "HP:0500066", "HP:0011003", "HP:0031730"]

    assert pheno_name_to_id(pheno_name) == expect


def test_pheno_id_to_term_HP():
    pheno_id = "HP:0000545"
    expect = "Myopia"

    assert pheno_id_to_term(pheno_id) == expect


def test_pheno_id_to_term_Orphanet():
    pass


def test_pheno_id_to_term_EFO():
    pass


def test_disease_id_to_name():
    pass
