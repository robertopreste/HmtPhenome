#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
# from quart_openapi import Resource
# from flask_restplus import Resource
# from app.api.views import api
# from app.site.models import User
# from app.api.models import user_schema, user_schema_many

# ns = api.namespace("users", description="Retrieve data from the Users table.")


# @ns.deprecated
# @ns.route("/")
# class UserList(Resource):
#     def get(self):
#         """
#         Get all the entries in Users table.
#         It is not recommended to run this query, as it may load a huge number of entries and consequently slow down your browser. Will return a list of entries.
#         """
#         q = User.query.all()
#         return user_schema_many.jsonify(q)
#
#
# @ns.route("/<int:user_id>")
# class UserId(Resource):
#     def get(self, user_id):
#         """
#         Find the information associated with the specified user id.
#         Will return a single entry.
#         """
#         q = User.query.filter(User.id == user_id).first()
#         return user_schema.jsonify(q)
