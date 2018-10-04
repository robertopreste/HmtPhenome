#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os

basedir = os.path.abspath(os.path.dirname(__file__))

WTF_CSRF_ENABLED = True
SECRET_KEY = "secret_key"

SQLALCHEMY_TRACK_MODIFICATIONS = False

# database
SQLALCHEMY_DATABASE_URI = "sqlite:///" + os.path.join(basedir, "database.db")
SQLALCHEMY_MIGRATE_REPO = os.path.join(basedir, "db_repo")

# DA CONFIGURARE
# mail server
MAIL_SERVER = "mail.server.com"
MAIL_PORT = 25
MAIL_USE_TLS = True
MAIL_USERNAME = "user.name"
MAIL_PASSWORD = "password1234"

# administrators
ADMINS = ["admin@mail.com"]
