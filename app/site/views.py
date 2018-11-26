#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
# import requests
# import json
from quart import Blueprint, render_template, request, redirect, url_for
from app.static import dbdata
from app.site.forms import QueryVariantsForm, QueryGenesForm, QueryPhenosForm, QueryDiseasesForm
from app.site.scripts import get_gene_from_variant, get_pheno_from_variant, get_diseases_from_variant, get_vars_from_gene_name, get_diseases_from_gene_name, get_genes_from_phenotype, get_vars_from_phenotype, get_diseases_from_phenotype, get_genes_from_disease_name, disease_id_to_name, get_vars_from_disease_name
# from flask import Blueprint, render_template, flash, redirect, session, url_for, request, g, jsonify, send_file
# from werkzeug.urls import url_parse

www = Blueprint("site", __name__)

# from sqlalchemy import or_, and_
# from app import app, db
# from config import ADMINS
# from .forms import LoginForm, RegistrationForm
# from .models import User


# Home Page
@www.route("/index", methods=["GET"])
@www.route("/home", methods=["GET"])
@www.route("/", methods=["GET"])
async def index():
    return await render_template("index.html",
                                 title="Home",
                                 latest_update=dbdata.latest_update)


@www.route("/about", methods=["GET"])
async def about():
    return await render_template("about.html",
                                 title="About")


@www.route("/contacts", methods=["GET"])
async def contacts():
    return await render_template("contacts.html",
                                 title="Contacts")


@www.route("/query", methods=["GET", "POST"])
async def query():
    if request.method == "GET":
        return await render_template("query.html",
                                     title="Query")

    elif request.method == "POST":
        form_var = QueryVariantsForm()
        form_gene = QueryGenesForm()
        form_phen = QueryPhenosForm()
        form_dis = QueryDiseasesForm()

        return redirect(url_for("site.results",
                                variant_input=form_var.variant_input.data,
                                variant_submit=form_var.variant_submit.data,
                                gene_chr=form_gene.gene_chr.data,
                                gene_input=form_gene.gene_input.data,
                                gene_submit=form_gene.gene_submit.data,
                                pheno_input=form_phen.pheno_input.data,
                                pheno_submit=form_phen.pheno_submit.data,
                                disease_input=form_dis.disease_input.data,
                                disease_submit=form_dis.disease_submit.data))


@www.route("/results", methods=["GET"])
async def results():

    variant_input = request.args.get("variant_input", "", type=str)
    variant_submit = request.args.get("variant_submit")
    gene_chr = request.args.get("gene_chr", "", type=str)
    gene_input = request.args.get("gene_input", "", type=str)
    gene_submit = request.args.get("gene_submit")
    pheno_input = request.args.get("pheno_input", "", type=str)
    pheno_submit = request.args.get("pheno_submit")
    disease_input = request.args.get("disease_input", "", type=str)
    disease_submit = request.args.get("disease_submit")

    if variant_submit == "True":
        var_chrom, var_rest = variant_input.split(":")
        if "-" in var_rest:
            var_start, var_end = var_rest.split("-")
        else:
            var_end = var_start = var_rest

        genes_df = get_gene_from_variant(var_chrom, var_start, var_end)
        phenos_df = get_pheno_from_variant(var_chrom, var_start, var_end)
        disease_df = get_diseases_from_variant(var_chrom, var_start, var_end)

        # TODO: move the following part to scripts, creating a proper class
        full_df = (disease_df.set_index("disease")
                   .join(phenos_df.set_index("phenotype"))
                   .reset_index())
        full_df["variant"] = full_df["ref_allele"] + full_df["start_pos"].astype(str) + full_df["alt_allele"]
        full_df = full_df[["ensembl_gene_id", "gene_name", "variation", "variant", "disease",
                           "phenotypes"]]
        # TODO: create a `phenotype_names` column based on data from `phenotypes`

    elif gene_submit == "True":
        vars_df = get_vars_from_gene_name(gene_input)
        phenos_df = get_diseases_from_gene_name(gene_input)  # TODO: get only phenotypes
        disease_df = get_diseases_from_gene_name(gene_input)  # TODO: get only diseases

    elif pheno_submit == "True":
        genes_df = get_genes_from_phenotype(pheno_input)  # TODO: get only genes
        vars_df = get_vars_from_phenotype(pheno_input)
        disease_df = get_diseases_from_phenotype(pheno_input)

    elif disease_submit == "True":
        disease_input_name = disease_id_to_name(disease_input)
        genes_df = get_genes_from_disease_name(disease_input_name)
        vars_df = get_vars_from_disease_name(disease_input_name)  # TODO: restrict results to selected disease
        # phenos_df = # TODO

    return await render_template("results.html",
                                 title="Results", full_df=full_df)


@www.errorhandler(404)
async def page_not_found(e):
    return await render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
async def internal_server_error(e):
    return await render_template("500.html", title="Error 500"), 500

