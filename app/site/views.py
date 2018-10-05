#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import requests
import json
from flask import Blueprint, render_template, flash, redirect, session, url_for, request, g, jsonify, send_file
from werkzeug.urls import url_parse

www = Blueprint("site", __name__)

from sqlalchemy import or_, and_
from app import app, db
from config import ADMINS
from .forms import LoginForm, RegistrationForm
from .models import User
from app.static import dbdata


# Home Page
@www.route("/index", methods=["GET"])
@www.route("/home", methods=["GET"])
@www.route("/", methods=["GET"])
def index():
    return render_template("index.html",
                           title="HmtPhenome",
                           latest_update=dbdata.latest_update)


@www.errorhandler(404)
def page_not_found(e):
    return render_template("404.html", title="Error 404"), 404


@www.errorhandler(500)
def internal_server_error(e):
    return render_template("500.html", title="Error 500"), 500

