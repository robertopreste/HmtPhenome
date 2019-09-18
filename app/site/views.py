#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import pprint
import os
import async_timeout
import pandas as pd
from quart import Blueprint, render_template, request, redirect, url_for, flash, session, send_from_directory
from app.static import dbdata
from app.site.forms import QueryVariantsForm, QueryGenesForm, QueryPhenosForm, \
    QueryDiseasesForm
from app.site.scripts import json_from_variant, network_from_variant_json, \
    json_from_gene, network_from_gene_json, json_from_phenotype, \
    network_from_phenotype_json, json_from_disease, network_from_disease_json, \
    parse_variant_string, parse_variant_HGVS, fallback_variant, create_dataframes

www = Blueprint("site", __name__)
dl_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "dls")


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
                                variant_pos=form_var.variant_pos.data,
                                variant_alt=form_var.variant_alt.data,
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
    variant_pos = request.args.get("variant_pos", "", type=str)
    variant_alt = request.args.get("variant_alt", "", type=str)
    variant_submit = request.args.get("variant_submit")
    gene_input = request.args.get("gene_input", "", type=str)
    gene_submit = request.args.get("gene_submit")
    pheno_input = request.args.get("pheno_input", "", type=str)
    pheno_submit = request.args.get("pheno_submit")
    disease_input = request.args.get("disease_input", "", type=str)
    disease_submit = request.args.get("disease_submit")

    if variant_submit == "True":
        res_type = "variant"
        # TODO: change this to the proper HGVS variant format

        variant_start = int(variant_pos)
        if variant_alt:
            res_el = ":".join([variant_chr, variant_pos + ">" + variant_alt.upper()])
            variant_end = int(variant_pos) + len(variant_alt) - 1
        else:
            res_el = ":".join([variant_chr, variant_pos])
            variant_end = variant_start
        # if "-" in variant_input:
        #     variant_start, variant_end = variant_input.split("-")
        # else:
        #     variant_end = variant_start = variant_input

        try:
            if variant_alt:
                variant_alt = variant_alt.upper()
            async with async_timeout.timeout(60) as cm:
                json_data = await json_from_variant(variant_chr, variant_start,
                                                    variant_alt, variant_end)
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
            async with async_timeout.timeout(120) as cm:
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
                                 "diseases": [], "phenotypes": []}
                    await flash("No results found!")
        except:
            json_data = {"variants": [], "genes": [],
                         "diseases": [], "phenotypes": []}
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
                                 "diseases": [], "phenotypes": []}
                    await flash("No results found!")
        except:
            json_data = {"variants": [], "genes": [],
                         "diseases": [], "phenotypes": []}
            await flash("No results found!")
        networks = network_from_disease_json(json_data)

    session["res_type"] = res_type
    session["res_el"] = res_el
    session["nodes"] = networks["nodes"]
    session["edges"] = networks["edges"]

    try:
        os.remove("app/site/dls/hmtphenome_data.xlsx")
    except:
        pass

    dfs = create_dataframes(json_data)
    with pd.ExcelWriter("app/site/dls/hmtphenome_data.xlsx") as writer:
        dfs["variants"].to_excel(writer, sheet_name="variants")
        dfs["genes"].to_excel(writer, sheet_name="genes")
        dfs["diseases"].to_excel(writer, sheet_name="diseases")
        dfs["phenotypes"].to_excel(writer, sheet_name="phenotypes")

    return await render_template("results.html",
                                 title="Results",
                                 res_type=res_type,
                                 res_el=res_el,
                                 # json_data=json_data,
                                 json_pretty=pprint.pformat(json_data),
                                 nodes=networks["nodes"],
                                 edges=networks["edges"],
                                 variants=networks["variants"],
                                 genes=networks["genes"],
                                 diseases=networks["diseases"],
                                 phenotypes=networks["phenotypes"],
                                 # parse_variant_string=parse_variant_string)
                                 parse_variant_string=parse_variant_HGVS)

    # elif request.method == "POST":
    #     return redirect(url_for("site.network"))


@www.route("/download_data", methods=["GET"])
async def download_data():
    # json_data = request.args.get("json_data", "")
    # dfs = create_dataframes(json_data)
    # with pd.ExcelWriter("app/site/dls/hmtphenome_data.xlsx") as writer:
    #     dfs["variants"].to_csv(writer, sheet_name="variants")
    #     dfs["genes"].to_csv(writer, sheet_name="genes")
    #     dfs["diseases"].to_csv(writer, sheet_name="diseases")
    #     dfs["phenotypes"].to_csv(writer, sheet_name="phenotypes")

    # @after_this_request
    # async def remove_file(response):
    #     os.remove("app/site/dls/hmtphenome_data.xlsx")
    #     return await response

    return await send_from_directory(dl_dir, file_name="hmtphenome_data.xlsx")

    # return await send_file("app/site/dls/hmtphenome_data.xlsx", mimetype="text/plain",
    #                        as_attachment=True,
    #                        attachment_filename="hmtphenome_data.xlsx")


@www.route("/network", methods=["GET"])
async def network():
    res_type = session.get("res_type")
    res_el = session.get("res_el")
    nodes = session.get("nodes")
    edges = session.get("edges")

    return await render_template("network.html",
                                 title="Network View",
                                 res_type=res_type,
                                 res_el=res_el,
                                 nodes=nodes,
                                 edges=edges)


@www.errorhandler(404)
async def page_not_found(e):
    return await render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
async def internal_server_error(e):
    return await render_template("500.html", title="Error 500"), 500

