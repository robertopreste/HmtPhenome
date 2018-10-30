#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
# from datetime import datetime
# from threading import Thread
# from app import db
# from app import db, login, mail, app
# from flask import render_template
# from flask_login import UserMixin
# from flask_mail import Message
# from werkzeug.security import generate_password_hash, check_password_hash
# from config import ADMINS


# class Diseases(db.Model):
#     __tablename__ = "Diseases"
#
#     id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
#     disease_sym = db.Column(db.String, index=True, unique=True, nullable=False)
#     disease_name = db.Column(db.String, index=True, unique=True, nullable=False)
#
#     def __repr__(self):
#         return """Diseases(id: {self.id}, disease_sym: {self.disease_sym}, disease_name: {self.disease_name})""".format(self=self)


# class Mitocarta(db.Model):
#     __tablename__ = "Mitocarta"
#
#     id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
#     gene_id = db.Column(db.Integer, index=True, unique=True, nullable=False)
#     ensembl_id = db.Column(db.String, nullable=True, default=None)
#     gene_symbol = db.Column(db.String, index=True, unique=True, nullable=False)
#     description = db.Column(db.Text, nullable=False)
#     hg_chr = db.Column(db.String, nullable=False)
#     hg_start = db.Column(db.Integer, nullable=False)
#     hg_end = db.Column(db.Integer, nullable=False)
#
#     def __repr__(self):
#         return """Genes(id: {self.id}, gene_id: {self.gene_id}, ensembl_id: {self.ensembl_id}, gene_symbol: {self.gene_symbol}, description: {self.description}, hg_chr: {self.hg_chr}, hg_start: {self.hg_start}, hg_end: {self.hg_end})""".format(self=self)


# class Phenotypes(db.Model):
#     __tablename__ = "Phenotypes"
#
#     id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
#     hpo_id = db.Column(db.String, index=True, unique=True, nullable=False)
#     hpo_term = db.Column(db.String, index=True, unique=True, nullable=False)
#
#     def __repr__(self):
#         return """Phenotypes(id: {self.id}, hpo_id: {self.hpo_id}, hpo_term: {self.hpo_term})""".format(self=self)


# class Variants(db.Model):
#     __tablename__ = "Variants"
#
#     id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
#     variant_id = db.Column(db.Integer, nullable=False)
#     group = db.Column(db.String, nullable=False)
#     locus = db.Column(db.String, index=True, unique=False, nullable=False)
#     variant = db.Column(db.String, index=True, unique=False, nullable=False)
#     nt_start = db.Column(db.Integer, nullable=False)
#     ref_rCRS = db.Column(db.String, nullable=False)
#     alt = db.Column(db.String, nullable=False)
#     nt_end = db.Column(db.Integer, nullable=False)
#     aa_change = db.Column(db.String, nullable=True, default=None)
#     codon_position = db.Column(db.Integer, nullable=True, default=None)
#     pathogenicity = db.Column(db.String, nullable=True, default=None)
#     disease_score = db.Column(db.Float, nullable=True, default=None)
#     clinvar = db.Column(db.String, nullable=True, default=None)
#     omim_num = db.Column(db.Integer, nullable=True, default=None)
#     omim_str = db.Column(db.String, nullable=True, default=None)
#     dbSNP = db.Column(db.String, nullable=True, default=None)
#     phenos = db.Column(db.String, nullable=True, default=None)
#
#     def __repr__(self):
#         return """Variants(id: {self.id}, variant_id: {self.variant_id}, group: {self.group}, locus: {self.locus}, variant: {self.variant}, nt_start: {self.nt_start}, ref_rCRS: {self.ref_rCRS}, alt: {self.alt}, nt_end: {self.nt_end}, aa_change: {self.aa_change}, codon_position: {self.codon_position}, pathogenicity: {self.pathogenicity}, disease_score: {self.disease_score}, clinvar: {self.clinvar}, omim_num: {self.omim_num}, omim_str: {self.omim_str}, dbNSP: {self.dbSNP}, phenos: {self.phenos})""".format(self=self)


# class Disease_vs_Pheno(db.Model):
#     __tablename__ = "Disease_vs_Pheno"
#
#     id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
#     disease_id = db.Column(db.String, index=True, nullable=False)
#     hpo_id = db.Column(db.String, index=True, nullable=False)
#
#     def __repr__(self):
#         return """Disease_vs_Pheno(id: {self.id}, disease_id: {self.disease_id}, hpo_id: {self.hpo_id})""".format(self=self)


# class User(UserMixin, db.Model):
#     __tablename__ = "User"
#
#     id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
#     username = db.Column(db.String, index=True, unique=True, nullable=False)
#     email = db.Column(db.String, index=True, unique=True, nullable=False)
#     first_name = db.Column(db.String, nullable=True, default=None)
#     last_name = db.Column(db.String, nullable=True, default=None)
#     password_hash = db.Column(db.String)
#     registr_date = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
#     last_access = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
#     approved = db.Column(db.Boolean, default=False)
#     downloads = db.relationship("Downloads", backref="User", lazy="dynamic")
#
#     def set_password(self, password):
#         self.password_hash = generate_password_hash(password)
#
#     def check_password(self, password):
#         return check_password_hash(self.password_hash, password)
#
#     def set_approval(self):
#         self.approved = True
#         db.session.commit()
#
#     def unset_approval(self):
#         self.approved = False
#         db.session.commit()
#
#     def update_last_access(self):
#         self.last_access = datetime.utcnow()
#         db.session.commit()
#
#     def __repr__(self):
#         return """User(id: {self.id}, username: {self.username})""".format(self=self)


# class Downloads(db.Model):
#     __tablename__ = "Downloads"
#
#     id = db.Column(db.Integer, nullable=False, autoincrement=True, primary_key=True)
#     dataset = db.Column(db.String, nullable=False)
#     dl_date = db.Column(db.DateTime, nullable=False, default=datetime.utcnow)
#     user_id = db.Column(db.Integer, db.ForeignKey("User.id"), nullable=False)
#
#     def __repr__(self):
#         return """Downloads(id: {self.id}, dataset: {self.dataset}, dl_date: {self.dl_date}, user_id: {self.user_id})""".format(self=self)


# @login.user_loader
# def load_user(id):
#     return User.query.get(int(id))
