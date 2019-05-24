#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import quart.flask_patch
import click
import datetime
import imp
import os
import pandas as pd
from quart import Quart
# from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_bootstrap import Bootstrap
from flask_marshmallow import Marshmallow
# from flask_cors import CORS
# from flask_login import LoginManager
# from flask_mail import Mail
from migrate.versioning import api
from config import SQLALCHEMY_DATABASE_URI
from config import SQLALCHEMY_MIGRATE_REPO


app = Quart(__name__)
# app = Flask(__name__)
# CORS(app)
app.config.from_object("config")
Bootstrap(app)
db = SQLAlchemy(app)
ma = Marshmallow(app)
# login = LoginManager(app)
# login.login_view = "site.login"
# mail = Mail(app)

# from app.api.views import res
from app.site.views import www

# from app.api.views import api
# from app.api.endpoints.users import ns as user_namespace


# api.add_namespace(user_namespace)

app.register_blueprint(www)
# app.register_blueprint(res)
# app.register_blueprint(www, static_folder="site/static")
# app.register_blueprint(res, url_prefix="/api")

from app.site.models import Mitocarta, Phenotypes


@app.cli.command()
def create_db():
    click.echo("Creating new database... ", nl=False)
    db.drop_all()
    db.create_all()
    click.echo("Done.")


@app.cli.command()
def migrate_db():
    from app.site.scripts import populate_genes, populate_phenos, \
        populate_diseases, populate_genes_autocomplete
    click.echo("Migrating database to new version... ", nl=False)

    with open("app/static/js/script.js", "w") as f:
        # f.write(populate_genes())
        f.write(populate_genes_autocomplete())
        f.write("\n")
        f.write(populate_phenos())
        f.write("\n")
        f.write(populate_diseases())
        f.write("\n")

    with open("app/static/dbdata.py", "w") as f:
        latest_update = "{} {}".format(datetime.date.today().strftime("%B"),
                                       str(datetime.date.today().year))
        f.write("latest_update = " + repr(latest_update) + "\n")

    click.echo("Done.")


@app.cli.command()
def update_db():
    click.echo("Updating database tables...")
    sources = ("Mitocarta", "HpoDisGenePhen", "HpoGenePhen", "HpoPhenGene",
               "Omim", "Orphanet", "GeneDiseaseAss", "VarDiseaseAss",
               "DiseaseMappings", "Phenotypes", "Diseases")
    for el in sources:
        click.echo("\tUpdating {} table... ".format(el), nl=False)
        df = pd.read_csv("app/update/data/tables/{}.csv".format(el))
        df.reset_index(drop=True, inplace=True)
        df.to_sql(name=el, con=db.engine, index=False, if_exists="replace")
        click.echo("Complete.")
    click.echo("Database correctly updated. Please migrate it before use.")



