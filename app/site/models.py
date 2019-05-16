#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from app import db


class Diseases(db.Model):
    __tablename__ = "Diseases"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    disease_id = db.Column(db.String(32), index=True, unique=True,
                           nullable=False)
    disease_name = db.Column(db.String(200), index=True, nullable=False)

    def __repr__(self):
        return """Diseases(id: {self.id}, disease_id: {self.disease_id}, 
        disease_name: {self.disease_name})""".format(self=self)


class Mitocarta(db.Model):
    __tablename__ = "Mitocarta"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    ensembl_id = db.Column(db.String(32), nullable=True, default=None)
    gene_symbol = db.Column(db.String(32), index=True, unique=True,
                            nullable=False)
    description = db.Column(db.Text, nullable=False)
    hg_chr = db.Column(db.String(32), nullable=False)
    hg_start = db.Column(db.Integer, nullable=False)
    hg_stop = db.Column(db.Integer, nullable=False)

    def __repr__(self):
        return """Mitocarta(id: {self.id}, ensembl_id: {self.ensembl_id}, 
        gene_symbol: {self.gene_symbol}, description: {self.description}, 
        hg_chr: {self.hg_chr}, hg_start: {self.hg_start}, 
        hg_stop: {self.hg_stop})""".format(self=self)


class Phenotypes(db.Model):
    __tablename__ = "Phenotypes"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    hpo_id = db.Column(db.String(32), index=True, unique=True, nullable=False)
    hpo_term_name = db.Column(db.String(200), index=True, unique=True,
                              nullable=False)

    def __repr__(self):
        return """Phenotypes(id: {self.id}, hpo_id: {self.hpo_id}, 
        hpo_term_name: {self.hpo_term_name})""".format(self=self)


class Omim(db.Model):
    __tablename__ = "Omim"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    mim_number = db.Column(db.Integer, nullable=False, index=True)
    mim_name = db.Column(db.String(200), nullable=False, index=True)
    prefix = db.Column(db.String(32), nullable=True, default=None)

    def __repr__(self):
        return """Omim(id: {self.id}, mim_number: {self.mim_number}, 
        mim_name: {self.mim_name}, prefix: {self.prefix})""".format(self=self)


class Orphanet(db.Model):
    __tablename__ = "Orphanet"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    orpha_num = db.Column(db.Integer, nullable=False, index=True)
    orpha_name = db.Column(db.String(200), nullable=False, index=True)

    def __repr__(self):
        return """Orphanet(id: {self.id}, orpha_num: {self.orpha_num}, 
        orpha_name: {self.orpha_name})""".format(self=self)


class GeneDiseaseAss(db.Model):
    __tablename__ = "GeneDiseaseAss"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    entrez_gene_id = db.Column(db.String(32), nullable=True, default=None)
    gene_symbol = db.Column(db.String(32), nullable=False, index=True)
    umls_disease_id = db.Column(db.String(32), nullable=False)
    disease_name = db.Column(db.String(200), nullable=False, index=True)
    score = db.Column(db.Float, nullable=True, default=None)

    def __repr__(self):
        return """GeneDiseaseAss(id: {self.id}, entrez_gene_id: {self.entrez_gene_id}, 
        gene_symbol: {self.gene_symbol}, umls_disease_id: {self.umls_disease_id}, 
        disease_name: {self.disease_name}, score: {self.score})""".format(self=self)


class VarDiseaseAss(db.Model):
    __tablename__ = "VarDiseaseAss"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    dbsnp_id = db.Column(db.String(32), nullable=False, index=True)
    umls_disease_id = db.Column(db.String(32), nullable=False)
    disease_name = db.Column(db.String(200), nullable=False, index=True)
    score = db.Column(db.Float, nullable=True, default=None)

    def __repr__(self):
        return """VarDiseaseAss(id: {self.id}, dbsnp_id: {self.dbsnp_id}, 
        umls_disease_id: {self.umls_disease_id}, disease_name: {self.disease_name}, 
        score: {self.score})""".format(self=self)


class DiseaseMappings(db.Model):
    __tablename__ = "DiseaseMappings"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    umls_disease_id = db.Column(db.String(32), nullable=False)
    disease_name = db.Column(db.String(200), nullable=False, index=True)
    vocabulary = db.Column(db.String(32), nullable=True, default=None)
    disease_id = db.Column(db.String(32), nullable=False, index=True)
    alt_disease_name = db.Column(db.String(400), nullable=False, index=True)

    def __repr__(self):
        return """DiseaseMappings(id: {self.id}, 
        umls_disease_id: {self.umls_disease_id}, 
        disease_name: {self.disease_name}, vocabulary: {self.vocabulary}, 
        disease_id: {self.disease_id}, 
        alt_disease_name: {self.alt_disease_name})""".format(self=self)


class HpoDisGenePhen(db.Model):
    __tablename__ = "HpoDisGenePhen"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    disease_id = db.Column(db.String(32), nullable=False, index=True)
    gene_symbol = db.Column(db.String(32), nullable=False, index=True)
    entrez_gene_id = db.Column(db.String(32), nullable=False)
    hpo_id = db.Column(db.String(32), nullable=False, index=True)
    hpo_term_name = db.Column(db.String(200), nullable=False, index=True)

    def __repr__(self):
        return """HpoDisGenePhen(id: {self.id}, disease_id: {self.disease_id}, 
        gene_symbol: {self.gene_symbol}, entrez_gene_id: {self.entrez_gene_id}, 
        hpo_id: {self.hpo_id}, hpo_term_name: {self.hpo_term_name})""".format(self=self)


class HpoGenePhen(db.Model):
    __tablename__ = "HpoGenePhen"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    entrez_gene_id = db.Column(db.String(32), nullable=False)
    gene_symbol = db.Column(db.String(32), nullable=False, index=True)
    hpo_term_name = db.Column(db.String(200), nullable=False, index=True)
    hpo_id = db.Column(db.String(32), nullable=False, index=True)

    def __repr__(self):
        return """HpoGenePhen(id: {self.id}, gene_id: {self.gene_id}, 
        gene_symbol: {self.gene_symbol}, hpo_term_name: {self.hpo_term_name}, 
        hpo_id: {self.hpo_id})""".format(self=self)


class HpoPhenGene(db.Model):
    __tablename__ = "HpoPhenGene"

    id = db.Column(db.Integer, nullable=False, autoincrement=True,
                   primary_key=True)
    hpo_id = db.Column(db.String(32), nullable=False, index=True)
    hpo_term_name = db.Column(db.String(200), nullable=False, index=True)
    entrez_gene_id = db.Column(db.String(32), nullable=False)
    gene_symbol = db.Column(db.String(32), nullable=False, index=True)

    def __repr__(self):
        return """HpoPhenGene(id: {self.id}, hpo_id: {self.hpo_id}, 
        hpo_term_name: {self.hpo_term_name}, gene_id: {self.gene_id}, 
        gene_symbol: {self.gene_symbol})""".format(self=self)

