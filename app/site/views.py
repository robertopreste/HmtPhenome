#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
# import requests
# import json
from quart import Blueprint, render_template
from app.static import dbdata
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
    return await render_template("query.html",
                                 title="Query")


@www.errorhandler(404)
async def page_not_found(e):
    return await render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
async def internal_server_error(e):
    return await render_template("500.html", title="Error 500"), 500

