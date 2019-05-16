#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
# import requests
# import json
import pprint
from quart import Blueprint, render_template, request, redirect, url_for, \
    jsonify, flash, escape
from app.static import dbdata
from app.site.forms import QueryVariantsForm, QueryGenesForm, QueryPhenosForm, \
    QueryDiseasesForm
from app.site.scripts import json_from_variant, network_from_variant_json, \
    json_from_gene, network_from_gene_json, json_from_phenotype, \
    network_from_phenotype_json, json_from_disease, network_from_disease_json


www = Blueprint("site", __name__)


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
                                 title="About HmtPhenome")


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
                                # gene_chr=form_gene.gene_chr.data,
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
    # gene_chr = request.args.get("gene_chr", "", type=str)
    gene_input = request.args.get("gene_input", "", type=str)
    gene_submit = request.args.get("gene_submit")
    pheno_input = request.args.get("pheno_input", "", type=str)
    pheno_submit = request.args.get("pheno_submit")
    disease_input = request.args.get("disease_input", "", type=str)
    disease_submit = request.args.get("disease_submit")

    if variant_submit == "True":
        if "-" in variant_input:
            variant_start, variant_end = variant_input.split("-")
        else:
            variant_end = variant_start = variant_input

        json_data = json_from_variant(variant_chr, variant_start, variant_end)
        networks = network_from_variant_json(json_data)
        if len(json_data["variants"]) == 0 and len(json_data["genes"]) == 0 \
                and len(json_data["diseases"]) == 0 and len(json_data["phenotypes"]) == 0:
            json_data = "{}"
            await flash("No results found!")

    elif gene_submit == "True":
        json_data = json_from_gene(gene_input)
        networks = network_from_gene_json(json_data)
        if len(json_data["variants"]) == 0 and len(json_data["genes"]) == 0 \
                and len(json_data["diseases"]) == 0 and len(json_data["phenotypes"]) == 0:
            json_data = "{}"
            await flash("No results found!")

    elif pheno_submit == "True":
        json_data = json_from_phenotype(pheno_input)
        networks = network_from_phenotype_json(json_data)
        if len(json_data["variants"]) == 0 and len(json_data["genes"]) == 0 \
                and len(json_data["diseases"]) == 0 and len(json_data["phenotype"]) == 0:
            json_data = "{}"
            await flash("No results found!")

    elif disease_submit == "True":
        json_data = json_from_disease(disease_input)
        networks = network_from_disease_json(json_data)
        if len(json_data["variants"]) == 0 and len(json_data["genes"]) == 0 \
                and len(json_data["diseases"]) == 0 and len(json_data["phenotype"]) == 0:
            json_data = "{}"
            await flash("No results found!")

    return await render_template("results.html",
                                 title="Results",
                                 json_data=pprint.pformat(json_data),
                                 nodes=networks["nodes"],
                                 edges=networks["edges"],
                                 variants=networks["variants"],
                                 genes=networks["genes"],
                                 diseases=networks["diseases"],
                                 phenotypes=networks["phenotypes"])


@www.errorhandler(404)
async def page_not_found(e):
    return await render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
async def internal_server_error(e):
    return await render_template("500.html", title="Error 500"), 500

