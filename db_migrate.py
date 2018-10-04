#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import imp
import datetime
from migrate.versioning import api
from app import db
from config import SQLALCHEMY_DATABASE_URI
from config import SQLALCHEMY_MIGRATE_REPO


v = api.db_version(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
migration = SQLALCHEMY_MIGRATE_REPO + ("/versions/%03d_migration.py" % (v + 1))
tmp_module = imp.new_module("old_model")
old_model = api.create_model(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)

exec(old_model, tmp_module.__dict__)

script = api.make_update_script_for_model(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO, tmp_module.meta, db.metadata)
open(migration, "wt").write(script)
api.upgrade(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)

v = api.db_version(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)

# expl: create the dbdata.py file with data to populate dropdown menus and details
with open("app/static/dbdata.py", "w") as d:

    # expl: display last update of the db in home page
    latest_update = "%s %s" % (datetime.date.today().strftime("%B"), str(datetime.date.today().year))
    d.write("latest_update = " + repr(latest_update) + "\n")


print("New migration saved as " + migration)
print("Current database version: " + str(v))
