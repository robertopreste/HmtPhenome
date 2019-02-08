#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask_wtf import FlaskForm
from wtforms import StringField, SelectField, SubmitField


class QueryVariantsForm(FlaskForm):
    variant_chr = SelectField("variant_chr", default=None)
    variant_input = StringField("variant_input", default=None)
    variant_submit = SubmitField("variant_submit")


class QueryGenesForm(FlaskForm):
    gene_chr = SelectField("gene_chr", default=None)
    gene_input = SelectField("gene_input", default=None)
    gene_submit = SubmitField("gene_submit")


class QueryPhenosForm(FlaskForm):
    pheno_input = StringField("pheno_input", default=None)
    pheno_submit = SubmitField("pheno_submit")


class QueryDiseasesForm(FlaskForm):
    disease_input = StringField("disease_input", default=None)
    disease_submit = SubmitField("disease_submit")

