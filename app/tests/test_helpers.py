#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pytest
from app.site.scripts import create_variant_string


def test_mt_snp_variant():
    chrom = "M"
    nt_start = 3308
    ref_all = "G"
    alt_all = "A"
    expect = "chr{}:{}{}_{}".format(chrom, nt_start, ref_all, alt_all)

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

