#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pprint
import async_timeout
from quart import Blueprint, render_template, request, redirect, url_for, flash
from app.static import dbdata
from app.site.forms import QueryVariantsForm, QueryGenesForm, QueryPhenosForm, \
    QueryDiseasesForm
from app.site.scripts import json_from_variant, network_from_variant_json, \
    json_from_gene, network_from_gene_json, json_from_phenotype, \
    network_from_phenotype_json, json_from_disease, network_from_disease_json, \
    parse_variant_string, fallback_variant

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
    gene_input = request.args.get("gene_input", "", type=str)
    gene_submit = request.args.get("gene_submit")
    pheno_input = request.args.get("pheno_input", "", type=str)
    pheno_submit = request.args.get("pheno_submit")
    disease_input = request.args.get("disease_input", "", type=str)
    disease_submit = request.args.get("disease_submit")

    if variant_submit == "True":
        res_type = "position"
        res_el = ":".join([variant_chr, variant_input])
        if "-" in variant_input:
            variant_start, variant_end = variant_input.split("-")
        else:
            variant_end = variant_start = variant_input

        try:
            async with async_timeout.timeout(60) as cm:
                json_data = await json_from_variant(variant_chr, variant_start, variant_end)
                # networks = network_from_variant_json(json_data)
                if cm.expired or len(json_data["variants"]) == 0:
                    json_data = {"variants": {}, "genes": {},
                                 "diseases": {}, "phenotypes": {}}
                    await flash("No results found!")
        except:
            json_data = {"variants": {}, "genes": {},
                         "diseases": {}, "phenotypes": {}}
            await flash("No results found!")
        networks = network_from_variant_json(json_data)

    elif gene_submit == "True":
        res_type = "gene"
        res_el = gene_input
        try:
            async with async_timeout.timeout(60) as cm:
                json_data = await json_from_gene(gene_input)
                # networks = network_from_gene_json(json_data)
                if cm.expired or len(json_data["genes"]) == 0:
                    json_data = {"variants": {}, "genes": {},
                                 "diseases": {}, "phenotypes": {}}
                    await flash("No results found!")
        except:
            json_data = {"variants": {}, "genes": {},
                         "diseases": {}, "phenotypes": {}}
            await flash("No results found!")
        networks = network_from_gene_json(json_data)

    elif pheno_submit == "True":
        res_type = "phenotype"
        res_el = pheno_input
        try:
            async with async_timeout.timeout(60) as cm:
                json_data = await json_from_phenotype(pheno_input)
                # networks = network_from_phenotype_json(json_data)
                if cm.expired or len(json_data["phenotypes"]) == 0:
                    json_data = {"variants": [], "genes": [],
                                 "diseases": [], "phenotype": []}
                    await flash("No results found!")
        except:
            json_data = {"variants": [], "genes": [],
                         "diseases": [], "phenotype": []}
            await flash("No results found!")
        networks = network_from_phenotype_json(json_data)

    elif disease_submit == "True":
        res_type = "disease"
        res_el = disease_input
        try:
            async with async_timeout.timeout(60) as cm:
                json_data = await json_from_disease(disease_input)
                # networks = network_from_disease_json(json_data)
                if cm.expired or len(json_data["diseases"]) == 0:
                    json_data = {"variants": [], "genes": [],
                                 "diseases": [], "phenotype": []}
                    await flash("No results found!")
        except:
            json_data = {"variants": [], "genes": [],
                         "diseases": [], "phenotype": []}
            await flash("No results found!")
        networks = network_from_disease_json(json_data)

    return await render_template("results.html",
                                 title="Results",
                                 res_type=res_type,
                                 res_el=res_el,
                                 json_data=pprint.pformat(json_data),
                                 nodes=networks["nodes"],
                                 edges=networks["edges"],
                                 variants=networks["variants"],
                                 genes=networks["genes"],
                                 diseases=networks["diseases"],
                                 phenotypes=networks["phenotypes"],
                                 parse_variant_string=parse_variant_string)


@www.errorhandler(404)
async def page_not_found(e):
    return await render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
async def internal_server_error(e):
    return await render_template("500.html", title="Error 500"), 500

