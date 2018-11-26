#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pandas as pd
from app.site.scripts import get_gene_from_variant


correct = pd.DataFrame({"ensembl_gene_id": "ENSG00000198888", "gene_name": "MT-ND1"}, index=[0])


def test_1():
    assert correct.equals(get_gene_from_variant("MT", "3308"))


def test_2():
    assert correct.equals(get_gene_from_variant("mt", "3308"))


def test_3():
    assert correct.equals(get_gene_from_variant("M", "3308"))


def test_4():
    assert correct.equals(get_gene_from_variant("m", "3308"))


def test_5():
    assert correct.equals(get_gene_from_variant("MT", 3308))


def test_6():
    assert correct.equals(get_gene_from_variant("mt", 3308))


def test_7():
    assert correct.equals(get_gene_from_variant("M", 3308))


def test_8():
    assert correct.equals(get_gene_from_variant("m", 3308))


