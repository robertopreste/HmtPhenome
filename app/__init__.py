#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_bootstrap import Bootstrap
from flask_marshmallow import Marshmallow
from flask_cors import CORS
from flask_login import LoginManager
from flask_mail import Mail

app = Flask(__name__)
CORS(app)
app.config.from_object("config")
Bootstrap(app)
db = SQLAlchemy(app)
ma = Marshmallow(app)
login = LoginManager(app)
login.login_view = "site.login"
mail = Mail(app)

from app.api.views import res
from app.site.views import www

from app.api.views import api
from app.api.endpoints.users import ns as user_namespace


api.add_namespace(user_namespace)

app.register_blueprint(www, static_folder="site/static")
app.register_blueprint(res, url_prefix="/api")
