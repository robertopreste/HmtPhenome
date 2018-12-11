#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
# import requests
# import json
import pprint
from quart import Blueprint, render_template, request, redirect, url_for, jsonify
from app.static import dbdata
from app.site.forms import QueryVariantsForm, QueryGenesForm, QueryPhenosForm, QueryDiseasesForm
from app.site.scripts import get_gene_from_variant, get_pheno_from_variant, get_diseases_from_variant, get_vars_from_gene_name, get_diseases_from_gene_name, get_genes_from_phenotype, get_vars_from_phenotype, get_diseases_from_phenotype, get_genes_from_disease_name, disease_id_to_name, get_vars_from_disease_name, final_from_variant, final_from_gene_name, network_from_variant, network_from_gene_name, network_from_gene, json_from_variant, network_from_variant_json
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
                                variant_chr=form_var.variant_chr.data,
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

    variant_chr = request.args.get("variant_chr", "", type=str)
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
        if "-" in variant_input:
            var_start, var_end = variant_input.split("-")
        else:
            var_end = var_start = variant_input

        # TODO: genes_df is useless, everything is found through disease_df and pheno_df
        genes_df = get_gene_from_variant(variant_chr, var_start, var_end)
        pheno_df = get_pheno_from_variant(variant_chr, var_start, var_end)
        disease_df = get_diseases_from_variant(variant_chr, var_start, var_end)

        # final_df = final_from_variant(genes_df, pheno_df, disease_df)
        # networks = network_from_variant(final_df)
        json_data = json_from_variant(disease_df, pheno_df)
        networks = network_from_variant_json(json_data)

    elif gene_submit == "True":
        vars_df = get_vars_from_gene_name(gene_input)
        # phenos_df = get_diseases_from_gene_name(gene_input)  # this is useless
        disease_df = get_diseases_from_gene_name(gene_input, True)
        if disease_df.shape[0] == 0:
            disease_df = get_diseases_from_gene_name(gene_input)
            # TEST
        #     final_df = final_from_gene_name(vars_df, disease_df)
        # else:
        #     final_df = final_from_gene_name(vars_df, disease_df, True)
        # networks = network_from_gene_name(final_df)
        # END TEST
        networks = network_from_gene(gene_input, vars_df, disease_df)

    elif pheno_submit == "True":
        genes_df = get_genes_from_phenotype(pheno_input)  # TODO: useless since it gets called when creating vars_df
        vars_df = get_vars_from_phenotype(pheno_input)
        disease_df = get_diseases_from_phenotype(pheno_input)

    elif disease_submit == "True":
        disease_input_name = disease_id_to_name(disease_input)
        genes_df = get_genes_from_disease_name(disease_input_name)
        vars_df = get_vars_from_disease_name(disease_input_name)  # TODO: restrict results to selected disease
        # phenos_df = # TODO

    return await render_template("results.html",
                                 title="Results", # final_df=final_df,
                                 json_data=pprint.pformat(json_data),
                                 nodes=networks["nodes"], edges=networks["edges"])


@www.errorhandler(404)
async def page_not_found(e):
    return await render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
async def internal_server_error(e):
    return await render_template("500.html", title="Error 500"), 500

