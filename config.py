#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os
from dotenv import load_dotenv

load_dotenv()

basedir = os.path.abspath(os.path.dirname(__file__))

WTF_CSRF_ENABLED = True
SECRET_KEY = os.getenv("SECRET_KEY") or "secret_key"

SQLALCHEMY_TRACK_MODIFICATIONS = False

# database
SQLALCHEMY_DATABASE_URI = "sqlite:///" + os.path.join(basedir, "hmtphenome.db")
SQLALCHEMY_MIGRATE_REPO = os.path.join(basedir, "db_repo")

# DA CONFIGURARE
# mail server
MAIL_SERVER = os.getenv("MAIL_SERVER") or "mail.server.com"
MAIL_PORT = os.getenv("MAIL_PORT") or 25
MAIL_USE_TLS = os.getenv("MAIL_USE_TLS") or True
MAIL_USERNAME = os.getenv("MAIL_USERNAME") or "user.name"
MAIL_PASSWORD = os.getenv("MAIL_PASSWORD") or "password1234"

# administrators
ADMINS = os.getenv("ADMINS") or "admin@mail.com"
