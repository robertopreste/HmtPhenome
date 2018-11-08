#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
# import requests
# import json
from quart import Blueprint, render_template, request, redirect, url_for
from app.static import dbdata
from app.site.forms import QueryVariantsForm, QueryGenesForm, QueryPhenosForm, QueryDiseasesForm
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
                                gene_chr=form_gene.gene_chr.data,
                                gene_input=form_gene.gene_input.data,
                                pheno_input=form_phen.pheno_input.data,
                                disease_input=form_dis.disease_input.data))


@www.route("/results", methods=["GET"])
async def results():

    variant_input = request.args.get("variant_input", "", type=str)
    gene_chr = request.args.get("gene_chr", "", type=str)
    gene_input = request.args.get("gene_input", "", type=str)
    pheno_input = request.args.get("pheno_input", "", type=str)
    disease_input = request.args.get("disease_input", "", type=str)

    return await render_template("results.html",
                                 title="Results")


@www.errorhandler(404)
async def page_not_found(e):
    return await render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
async def internal_server_error(e):
    return await render_template("500.html", title="Error 500"), 500

