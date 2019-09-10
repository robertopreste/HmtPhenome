#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from app.site.models import Mitocarta, Phenotypes, Diseases, Omim, Orphanet, \
    GeneDiseaseAss, VarDiseaseAss, DiseaseMappings, HpoDisGenePhen
import pandas as pd
import numpy as np
import requests
import re
import apybiomart as apy
from fuzzywuzzy import fuzz
import json
from typing import List, Union, Optional
import asyncio
import aiohttp
import async_timeout
from app.site.classes import Node


def get_genes() -> dict:
    """Retrieve Mitocarta genes from the database and create a dictionary that
    will be used to populate the related dropdown menu in query page.

    :return: dict
    """
    chr_dict = {}
    chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
            "chrX", "chrM"]
    for chrom in chrs:
        lista = set()
        q = Mitocarta.query.filter(Mitocarta.hg_chr == chrom).all()
        for el in q:
            lista.add((el.gene_symbol, el.ensembl_id))
        chr_dict[chrom] = sorted(list(lista))

    return chr_dict


def get_genes_autocomplete() -> list:
    """Retrieve genes from the database and create a list that will be used
    to populate the related dropdown menu in the query page.

    :return: list
    """
    gene_list = []
    q = Mitocarta.query.all()
    for el in q:
        if el.ensembl_id is not None:
            gene_list.append([el.gene_symbol, el.ensembl_id])

    return gene_list


def get_phenos() -> list:
    """Retrieve phenotypes from the database and create a list that will be used
    to populate the related dropdown menu in the query page.

    :return: list
    """
    pheno_list = []
    q = Phenotypes.query.all()
    for el in q:
        pheno_list.append([el.hpo_term_name, el.hpo_id])

    return pheno_list


def get_diseases() -> list:
    """Retrieve diseases from the database and create a list that will be used
    to populate the related dropdown menu in the query page.

    :return: list
    """
    disease_list = []
    q = Diseases.query.all()
    for el in q:
        disease_list.append([el.disease_name, el.disease_id])

    return disease_list


def populate_phenos() -> str:
    """Add phenotypes data for the autocomplete function in the script.js file,
    during the update of the database.

    :return: str
    """
    phenos = get_phenos()
    base_string = """
var pheno_compl = document.getElementById("pheno_input"); 
new Awesomplete(pheno_compl, {list: %s});     
    """ % repr(phenos)

    return base_string


def populate_diseases() -> str:
    """Add diseases data for the autocomplete function in the script.js file,
    during the update of the database.

    :return: str
    """
    diseases = get_diseases()
    base_string = """
var disease_compl = document.getElementById("disease_input"); 
new Awesomplete(disease_compl, {list: %s});     
    """ % repr(diseases)

    return base_string


def populate_genes_autocomplete() -> str:
    """Add genes for the autocomplete function in the script.js file,
    during the update of the database.

    :return: str
    """
    genes = get_genes_autocomplete()
    base_string = """
var gene_compl = document.getElementById("gene_input");
new Awesomplete(gene_compl, {list: %s});
    """ % repr(genes)

    return base_string


def populate_genes() -> str:
    """Add genes for each chromosome in the script.js file, during the update
    of the database, in order to populate the related dropdown menu in the
    query page.

    :return: str
    """
    base_string = """
function populateGenes(s1, s2) {
    // Populate the genes dropdown menu based on selected chromosome
    s1 = document.getElementById(s1);
    s2 = document.getElementById(s2);
    var optionArray;

    s2.innerHTML = "--All Genes--";
    """

    chr_dict = get_genes()
    for n, chrom in enumerate(chr_dict):
        if n == 0:
            base_string += """
    if (s1.value == "%s") {
        optionArray = ["A|--All Genes--", """ % chrom
            for el in chr_dict[chrom]:
                base_string += """
        "%s|%s", """ % (el[1], el[0])
        else:
            base_string += """]; 
    } else if (s1.value == "%s") {
        optionArray = ["A|--All Genes--", """ % chrom

            for el in chr_dict[chrom]:
                base_string += """
        "%s|%s", """ % (el[1], el[0])

    base_string += """]; 
    }

    for (var option in optionArray) {
        var pair = optionArray[option].split("|");
        var newOption = document.createElement("option");

        newOption.value = pair[0];
        newOption.innerHTML = pair[1];
        s2.options.add(newOption);
    }
}
"""

    return base_string


async def get_json_request(url, headers=None):
    """Asynchronous request to a generic url, returning a JSON dictionary.

    :param url: url to request

    :param headers: dict of headers
        (default: {"Content-Type": "application/json"})

    :return:
    """
    if headers is None:
        headers = {"Content-Type": "application/json"}
    async with aiohttp.ClientSession() as session:
        async with session.get(url,
                               ssl=False,
                               headers=headers) as res:
            return await res.json()


def get_json_sync_request(url, headers=None):
    """Synchronous request to a generic url, returning a JSON dictionary.

    :param url: url to request

    :param headers: dict of headers
        (default: {"Content-Type": "application/json"})

    :return:
    """
    if headers is None:
        headers = {"Content-Type": "application/json"}
    res = requests.get(url, headers=headers)
    return res.json()


async def pheno_name_to_id(pheno_name: str) -> List[str]:
    """Retrieve the related HP id from a given phenotype name.

    :param str pheno_name: phenotype name

    :return: List[str] with the related id(s)
    """
    url = "https://hpo.jax.org/api/hpo/search?q={}".format(pheno_name)
    res = await get_json_request(url)
    ids = []

    if res["termsTotalCount"] == 0:
        return ids
    elif res["termsTotalCount"] == 1:
        ids.append(res["terms"][0]["id"])
    else:
        ids = [el["id"] for el in res["terms"]]

    return ids


async def pheno_id_to_term(pheno_id: str) -> str:
    """Retrieve the common phenotype name from a given ID. It can retrieve
    the data either from the local database (for HP IDs) or from the web
    (for EFO IDs).

    :param str pheno_id: HP or EFO ID

    :return: str with the related common name
    """
    pheno_name = ""
    if pheno_id.startswith("HP:"):
        try:
            pheno_name = Phenotypes.query.filter(
                Phenotypes.hpo_id == pheno_id
            ).first().hpo_term_name
        except AttributeError:
            return pheno_name
    elif pheno_id.startswith("Orphanet:"):
        orpha_id = pheno_id.strip("Orphanet:")
        pheno_name = Orphanet.query.filter(
            Orphanet.orpha_num == orpha_id
        ).first().orpha_name
    elif pheno_id.startswith("EFO:"):
        efo_id = pheno_id.strip("EFO:")
        base_url = "https://www.ebi.ac.uk/ols/api/ontologies/efo/terms?"
        iri_url = "iri=http://www.ebi.ac.uk/efo/EFO_{}".format(efo_id)
        res = await get_json_request(base_url + iri_url)
        pheno_name = res["_embedded"]["terms"][0]["label"]

    return pheno_name


def disease_id_to_name(disease_id: str) -> str:
    """Convert a given disease ID into its common name.

    :param str disease_id: query disease ID

    :return: str related disease name
    """
    try:
        q = Diseases.query.filter(Diseases.disease_id == disease_id).first()
        return q.disease_name
    except AttributeError:
        try:
            if disease_id.startswith("OMIM:"):
                q = Omim.query.filter(
                    Omim.mim_number == int(disease_id.strip("OMIM:"))
                ).first()
                return q.mim_name
            elif disease_id.startswith("ORPHA:"):
                q = Orphanet.query.filter(
                    Orphanet.orpha_num == int(disease_id.strip("ORPHA:"))
                ).first()
                return q.orpha_name
        except AttributeError:
            return ""


def disease_name_to_id(disease_name: str) -> str:
    """Convert a given disease name into the related Omim or Orphanet ID.

    :param str disease_name: query disease name

    :return: str related disease ID from OMIM or ORPHANET
    """
    # TODO: this should be looked for in the DiseaseMappings table 
    try:
        q = Diseases.query.filter(Diseases.disease_name == disease_name).first()
        return q.disease_id
    except AttributeError:
        try:
            q = Omim.query.filter(Omim.mim_name == disease_name).first()
            return "OMIM:{}".format(q.mim_number)  # TODO: check also mim symbol
        except AttributeError:
            try:
                q = Orphanet.query.filter(
                    Orphanet.orpha_name == disease_name
                ).first()
                return "ORPHA:{}".format(q.orpha_num)
            except AttributeError:
                return ""


def umls_to_disease_id(umls_id: str) -> str:
    """Retrieve the disease ID associated to the given UMLS ID.

    Priority is given to OMIM IDs, then to others.

    :param str umls_id: input UMLS ID

    :return: str the related disease ID
    """
    q = DiseaseMappings.query.filter(
        DiseaseMappings.umls_disease_id == umls_id).all()
    # first try with Omim
    omims = list(filter(lambda x: x.vocabulary == "OMIM", q))
    if len(omims) > 0:
        for el in omims:
            dis_q = Omim.query.filter(Omim.mim_number == el.disease_id,
                                      Omim.prefix == "Number Sign").first()
            if dis_q:
                return "OMIM:{}".format(dis_q.mim_number)
    # second try with Orphanet
    orphas = list(filter(lambda x: x.vocabulary == "ORDO", q))
    if len(orphas) > 0:
        for el in orphas:
            dis_q = Orphanet.query.filter(Orphanet.orpha_num == el.disease_id).first()
            if dis_q:
                return "ORPHA:{}".format(dis_q.orpha_num)
    # third try with DO
    dos = list(filter(lambda x: x.vocabulary == "DO", q))
    if len(dos) > 0:
        return "DO:{}".format(dos[0].disease_id)
    # last try with Mondo
    mondos = list(filter(lambda x: x.vocabulary == "MONDO", q))
    if len(mondos) > 0:
        return list(mondos)[0].disease_id

    if q:
        return q[0].vocabulary + ":" + q[0].disease_id
    else:
        return ""


def create_variant_string(chrom: Union[int, str],
                          nt_start: Union[int, str],
                          ref_all: str,
                          alt_all: str) -> str:
    """Create a string with the standard variant format:
    chrX:start_positionREF_ALL>ALT_ALL.

    :param Union[int, str] chrom: chromosome name (1:22, X, Y, M)

    :param Union[int, str] nt_start: start position of the variant

    :param str ref_all: reference allele

    :param str alt_all: alternate allele

    :return: str with variant formatted according to current standards
    """
    if type(ref_all) is float or type(alt_all) is float:
        return "chr_:_>_"
    base_str = "chr{}:{}{}"
    try:
        change = "{}>{}".format(ref_all.upper(), alt_all.upper())
        if alt_all == "-":  # deletion
            change = "del{}".format(ref_all.upper())
        elif ref_all == "-":  # insertion
            change = "ins{}".format(alt_all.upper())
    except TypeError:
        chrom = nt_start = "_"
        change = "_>_"

    return base_str.format(chrom, nt_start, change)


def parse_variant_string(variant: str) -> set:
    """Parse the variant string created by create_variant_string().
    The input variant string should be in the format chr:posA>T.

    :param str variant: input variant string

    :return: list(var_chr, var_pos, var_ref, var_alt)
    """
    rgx = re.compile(r"(chr.+):(\d+)(\w+)")
    return rgx.findall(variant)[0]


def parse_variant_HGVS(variant: str) -> set:
    """Parse the variant string created by variant_to_HGVS().

    :param variant: input variant string
    :return: list(var_chr, var_pos, var_ref, var_alt)
    """
    rgx = re.compile(r"(NC_.+):[gm]\.(\d+_*\d*)([\w>]+)")
    return rgx.findall(variant)[0]


def create_variant_HGVS(chrom: Union[int, str],
                        nt_start: Union[int, str],
                        ref_all: str,
                        alt_all: str) -> str:
    """Convert the given variant string to HGVS format.

    :param chrom: chromosome name (1:22, X, Y, M)
    :param nt_start: start position of the variant
    :param ref_all: reference allele
    :param alt_all: alternate allele
    :return: str variant in HGVS format
    """
    base_str = "{}:{}{}"

    chroms = {"1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
              "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
              "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
              "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
              "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
              "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
              "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
              "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
              "M": "NC_012920.1", "MT": "NC_012920.1"}

    ref = chroms[str(chrom)]
    pos = "m.{}" if str(chrom) in ["M", "MT"] else "g.{}"
    try:
        if alt_all == "-":  # deletion
            pos = pos.format("{}_{}".format(nt_start, int(nt_start) + len(ref_all)))
            change = "del"
        elif ref_all == "-":  # insertion
            pos = pos.format("{}_{}".format(nt_start, int(nt_start) + 1))
            change = "ins{}".format(alt_all.upper())
        else:  # SNP
            pos = pos.format(nt_start)
            change = "{}>{}".format(ref_all.upper(), alt_all.upper())
    except TypeError:
        pos = change = "_"

    return base_str.format(ref, pos, change)



async def ensembl_gene_id_to_entrez(ens_gene_id: str) -> pd.DataFrame:
    """Convert an Ensembl gene ID to its related Entrez gene ID, using Biomart.

    :param str ens_gene_id: query Ensembl gene ID

    :return: pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
        "entrez_gene_id"])
    """
    res = await apy.aquery(attributes=["ensembl_gene_id", "external_gene_name",
                                       "entrezgene_id"],
                           filters={"link_ensembl_gene_id": ens_gene_id},
                           dataset="hsapiens_gene_ensembl")
    res.rename({"Gene stable ID": "ensembl_gene_id", "Gene name": "gene_name",
                "NCBI gene ID": "entrez_gene_id"}, axis=1, inplace=True)

    return res


# FROM VARIANT POSITION #


async def get_dbsnp_from_variant(chrom: Union[int, str],
                           var_start: Union[int, str],
                           var_end: Optional[Union[int, str]] = None) -> pd.DataFrame:
    """Find the dbSNP ID related to a given variant.

    :param Union[int, str] chrom: chromosome name (1:22, X, Y, MT)

    :param Union[int, str] var_start: variant starting position

    :param Optional[Union[int, str]] var_end: variant ending position

    :return: pd.DataFrame(columns=["dbsnp_id", "variant"])
    """
    chrom = str(chrom).lstrip("chr").upper()
    var_start = str(var_start)
    if chrom == "M":
        chrom = "MT"
    if var_end is None:
        var_end = var_start
    else:
        var_end = str(var_end)

    res = await apy.aquery(attributes=["chr_name", "chrom_start",
                                       "consequence_allele_string", "refsnp_id"],
                           filters={"chr_name": chrom, "start": var_start,
                                    "end": var_end},
                           dataset="hsapiens_snp")
    res.rename({"Chromosome/scaffold name": "chromosome",
                "Chromosome/scaffold position start (bp)": "start_pos",
                "Consequence specific allele": "ref/alt allele",
                "Variant name": "dbsnp_id"}, axis=1, inplace=True)
    if res.shape[0] != 0:
        res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
        res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
        res = res[res["alt_allele"] != "HGMD_MUTATION"]
        res["variant"] = [create_variant_HGVS(el.chromosome, el.start_pos,
                                                el.ref_allele, el.alt_allele)
                          for el in res.itertuples()]
        res.drop_duplicates(subset="variant", inplace=True)
        res.drop(["chromosome", "start_pos", "ref_allele", "alt_allele",
                  "ref/alt allele"], axis=1, inplace=True)
    else:
        return pd.DataFrame(columns=["dbsnp_id", "variant"])

    return res


def get_gene_from_variant(chrom: Union[int, str],
                          var_start: Union[int, str],
                          var_end: Optional[Union[int, str]] = None) \
        -> pd.DataFrame:
    """Retrieve the gene to which the provided variant belongs.

    :param Union[int, str] chrom: chromosome name (chr + 1:22, X, Y, M)

    :param Union[int, str] var_start: variant starting position

    :param Optional[Union[int, str]] var_end: variant ending position

    :return: pd.DataFrame(columns=["ensembl_gene_id", "gene_name"])
    """
    if chrom == "MT":
        chrom = "M"
    if not chrom.startswith("chr"):
        chrom = "chr{}".format(chrom)
    if var_end is None:
        var_end = var_start

    q = Mitocarta.query.filter(Mitocarta.hg_chr == chrom,
                               Mitocarta.hg_start <= var_start,
                               Mitocarta.hg_stop >= var_end).first()
    try:
        res = pd.DataFrame({"ensembl_gene_id": [q.ensembl_id],
                            "gene_name": [q.gene_symbol]})
    except AttributeError:
        res = pd.DataFrame(columns=["ensembl_gene_id", "gene_name"])

    return res


def get_diseases_from_dbsnp(dbsnp_id: str) -> pd.DataFrame:
    """Retrieve diseases associated with the given dbSNP ID.

    :param str dbsnp_id: query dbSNP ID

    :return: pd.DataFrame(columns=["dbsnp_id", "umls_disease_id",
        "disease_id" "disease_name", "ass_score"])
    """
    q = VarDiseaseAss.query.filter(VarDiseaseAss.dbsnp_id == dbsnp_id).all()
    df = pd.DataFrame(columns=["dbsnp_id", "umls_disease_id", "disease_id",
                               "disease_name", "ass_score"])
    if q:
        for el in q:
            new_row = pd.DataFrame({"dbsnp_id": [el.dbsnp_id],
                                    "umls_disease_id": [el.umls_disease_id],
                                    "disease_id": [umls_to_disease_id(el.umls_disease_id)],
                                    "disease_name": [el.disease_name],
                                    "ass_score": [el.score]})
            df = df.append(new_row, ignore_index=True)
        df.sort_values(by="ass_score", ascending=False, inplace=True)

    return df


def get_phenos_from_umls(umls_id: str) -> pd.DataFrame:
    """Retrieve phenotypes associated with the given disease using its UMLS ID.

    :param str umls_id: query UMLS ID

    :return: pd.DataFrame(columns=["umls_disease_id", "disease_name",
        "disease_id", "phenotype_id", "phenotype_name"])
    """
    q = DiseaseMappings.query.filter(
        DiseaseMappings.umls_disease_id == umls_id,
        DiseaseMappings.vocabulary.in_(["OMIM", "ORPHA"])
    ).all()
    df = pd.DataFrame(columns=["umls_disease_id", "disease_name", "disease_id",
                               "phenotype_id", "phenotype_name"])
    if q:
        for el in q:
            disease_id = "{}:{}".format(el.vocabulary, el.disease_id)
            phenos = HpoDisGenePhen.query.filter(
                HpoDisGenePhen.disease_id == disease_id
            ).all()
            if phenos:
                for ph in phenos:
                    new_row = pd.DataFrame(
                        {"umls_disease_id": [el.umls_disease_id],
                         "disease_name": [el.disease_name],
                         "disease_id": [disease_id],
                         "phenotype_id": [ph.hpo_id],
                         "phenotype_name": [ph.hpo_term_name]}
                    )
                    df = df.append(new_row, ignore_index=True)

    df.drop_duplicates(inplace=True)

    return df


async def json_from_variant(variant_chr: Union[int, str],
                            variant_start: Union[int, str],
                            variant_end: Optional[Union[int, str]] = None) -> dict:
    """Create the final json structure from variant data.

    :param Union[int, str] variant_chr: chromosome name (chr + 1:22, X, Y, M)

    :param Union[int, str] variant_start: variant starting position

    :param Optional[Union[int, str]] variant_end: variant ending position

    :return: dict json("variants": [variants list])
    """
    gene = get_gene_from_variant(variant_chr, variant_start, variant_end)
    dbsnps = await get_dbsnp_from_variant(variant_chr, variant_start, variant_end)
    disease_df = pd.DataFrame(columns=["dbsnp_id", "umls_disease_id",
                                       "disease_id", "disease_name", "ass_score"])
    phenos_df = pd.DataFrame(columns=["umls_disease_id", "disease_name",
                                      "disease_id", "phenotype_id",
                                      "phenotype_name"])

    gene_json = json.loads(gene.to_json(orient="records"))

    if dbsnps.shape[0] != 0:
        dbsnp_id_unique = dbsnps.dbsnp_id.unique()
        for el in dbsnp_id_unique:
            disease_df = disease_df.append(get_diseases_from_dbsnp(el),
                                           ignore_index=True)

    disease_json = json.loads(disease_df.to_json(orient="records"))

    if disease_df.shape[0] != 0:
        umls_id_unique = disease_df.umls_disease_id.unique()
        for el in umls_id_unique:
            phenos_df = phenos_df.append(get_phenos_from_umls(el),
                                         ignore_index=True)

    phenos_json = json.loads(phenos_df.to_json(orient="records"))

    df = pd.DataFrame(columns=["variant", "dbsnp_id", "gene_name",
                               "ensembl_gene_id", "umls_disease_id",
                               "disease_id",
                               "disease_name", "phenotype_id",
                               "phenotype_name"])
    for var in dbsnps.itertuples():
        rel_dis = disease_df[disease_df.dbsnp_id == var.dbsnp_id]
        if rel_dis.shape[0] != 0:
            for dis in rel_dis.itertuples():
                rel_phenos = phenos_df[
                    phenos_df.umls_disease_id == dis.umls_disease_id
                    ]
                if rel_phenos.shape[0] != 0:
                    for phen in rel_phenos.itertuples():
                        new_row = pd.DataFrame(
                            {"variant": [var.variant],
                             "dbsnp_id": [var.dbsnp_id],
                             "gene_name": [gene.gene_name.get(0, "none")],
                             "ensembl_gene_id": [gene.ensembl_gene_id.get(
                                 0, "none"
                             )],
                             "umls_disease_id": [dis.umls_disease_id],
                             "disease_id": [dis.disease_id],
                             "disease_name": [dis.disease_name],
                             "phenotype_id": [phen.phenotype_id],
                             "phenotype_name": [phen.phenotype_name]})
                        df = df.append(new_row, ignore_index=True)
                else:
                    new_row = pd.DataFrame(
                        {"variant": [var.variant], "dbsnp_id": [var.dbsnp_id],
                         "gene_name": [gene.gene_name.get(0, "none")],
                         "ensembl_gene_id": [gene.ensembl_gene_id.get(0,
                                                                      "none")],
                         "umls_disease_id": [dis.umls_disease_id],
                         "disease_id": [dis.disease_id],
                         "disease_name": [dis.disease_name],
                         "phenotype_id": [""], "phenotype_name": [""]})
                    df = df.append(new_row, ignore_index=True)
        else:
            new_row = pd.DataFrame(
                {"variant": [var.variant], "dbsnp_id": [var.dbsnp_id],
                 "gene_name": [gene.gene_name.get(0, "none")],
                 "ensembl_gene_id": [gene.ensembl_gene_id.get(0, "none")],
                 "umls_disease_id": [""], "disease_id": [""], "disease_name": [""],
                 "phenotype_id": [""], "phenotype_name": [""]})
            df = df.append(new_row, ignore_index=True)

    vars_json = json.loads(df.to_json(orient="records"))

    final_json = dict()
    final_json["variants"] = vars_json
    final_json["genes"] = gene_json
    final_json["diseases"] = disease_json
    final_json["phenotypes"] = phenos_json

    return final_json


async def fallback_variant(variant_chr: Union[int, str],
                           variant_start: Union[int, str],
                           variant_end: Optional[Union[int, str]] = None) -> dict:
    """Fallback function, called when the final json structure from 
    variant data is empty.

    :param Union[int, str] variant_chr: chromosome name (chr + 1:22, X, Y, M)

    :param Union[int, str] variant_start: variant starting position

    :param Optional[Union[int, str]] variant_end: variant ending position

    :return: dict json("variants": [variants list])
    """
    dbsnps = await get_dbsnp_from_variant(variant_chr, variant_start, variant_end)
    final_json = {"variants": [{"variant": el.variant, "dbsnp_id": el.dbsnp_id}
                               for el in dbsnps.itertuples()],
                  "genes": [], "diseases": [], "phenotypes": []}
    return final_json


def network_from_variant_json(final_json: dict) -> dict:
    """Create the required nodes and edges dictionaries to build the network
    from variant data.

    :param dict final_json: output from json_from_variant()

    :return: dict("nodes": [nodes list], "edges": [edges list])
    """
    var_json = final_json["variants"]
    variants = []
    variants.extend([el["variant"] for el in var_json])
    variants = set(variants)
    variants.discard("")

    genes = []
    genes.extend([el["gene_name"] for el in var_json
                  if el["gene_name"] != "none"])
    genes = set(genes)
    genes.discard("")

    diseases = []
    diseases.extend([el["disease_name"] for el in var_json])
    diseases = set(diseases)
    diseases.discard("")

    phenotypes = []
    phenotypes.extend([el["phenotype_name"] for el in var_json])
    phenotypes = set(phenotypes)
    phenotypes.discard("")

    nodes = []
    edges = []
    ids = 0
    id_dict = {}
    for el in variants:
        ids += 1
        node = Node("v", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el, "color": {"background": "#F9CF45",
                                                        "border": "#CCAA39"}})
    for el in genes:
        ids += 1
        node = Node("g", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el, "color": {"background": "#739E82",
                                                        "border": "#5F826B"}})
    for el in diseases:
        if el != "":
            ids += 1
            node = Node("d", ids, el)
            id_dict[node.uid] = ids
            nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                            "border": "#B06A57"}})
    for el in phenotypes:
        if el != "":
            ids += 1
            node = Node("p", ids, el)
            id_dict[node.uid] = ids
            nodes.append({"id": ids, "label": el,
                          "color": {"background": "#93B5C6",
                                    "border": "#7995A3"}})

    connected_nodes = set()
    vars_set = set()  # (variant, dbsnp_id, gene_name)
    gene_set = set()  # (gene_name, ensembl_gene_id)
    dise_set = set()  # (disease_name, umls_disease_id)
    phen_set = set()  # (phenotype_name, phenotype_id)

    for el in var_json:
        # variant to gene
        try:
            if genes:  # avoid "none" genes
                edges.append({
                    "from": id_dict[hash("v" + el["variant"])],
                    "to": id_dict[hash("g" + el["gene_name"])]
                })
                connected_nodes.add(id_dict[hash("v" + el["variant"])])
                vars_set.add((el["variant"], el["dbsnp_id"], el["gene_name"]))
                connected_nodes.add(id_dict[hash("g" + el["gene_name"])])
                gene_set.add((el["gene_name"], el["ensembl_gene_id"]))
        except:
            pass
        # variant to diseases
        try:
            if diseases:
                edges.append({
                    "from": id_dict[hash("v" + el["variant"])],
                    "to": id_dict[hash("d" + el["disease_name"])]
                })
                connected_nodes.add(id_dict[hash("v" + el["variant"])])
                vars_set.add((el["variant"], el["dbsnp_id"], el["gene_name"]))
                connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
                dise_set.add((el["disease_name"], el["disease_id"]))
        except:
            pass
        # disease to phenotypes
        try:
            if el["phenotype_name"] != "" \
                    and id_dict[hash("d" + el["disease_name"])] != id_dict[hash("p" + el["phenotype_name"])]:
                edges.append({
                    "from": id_dict[hash("d" + el["disease_name"])],
                    "to": id_dict[hash("p" + el["phenotype_name"])]
                })
                connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
                dise_set.add((el["disease_name"], el["disease_id"]))
                connected_nodes.add(id_dict[hash("p" + el["phenotype_name"])])
                phen_set.add((el["phenotype_name"], el["phenotype_id"]))
        except:
            pass

    # delete orphan nodes
    true_nodes = [el for el in nodes if el["id"] in connected_nodes]

    return {"nodes": true_nodes, "edges": edges, "variants": vars_set,
            "genes": gene_set, "diseases": dise_set, "phenotypes": phen_set}


# FROM GENE #


async def get_vars_from_gene(ens_gene_id: str) -> pd.DataFrame:
    """Retrieve all variants associated with a specific Ensembl gene ID,
    using Ensembl REST.

    :param str ens_gene_id: query Ensembl gene ID

    :return: pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
        "chromosome", "ref_allele", "start_pos", "alt_allele", "variant",
        "dbsnp_id"])
    """
    try:
        gene_name = Mitocarta.query.filter(
            Mitocarta.ensembl_id == ens_gene_id
        ).first().gene_symbol
    except AttributeError:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                     "chromosome", "ref_allele", "start_pos",
                                     "alt_allele", "variant", "dbsnp_id"])

    server = "https://rest.ensembl.org"
    ext = "/overlap/id/{}?feature=variation".format(ens_gene_id)
    res = await get_json_request(server + ext)

    df = pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "chromosome",
                               "ref_allele", "start_pos", "alt_allele",
                               "variant", "dbsnp_id"])
    if len(res) == 0:
        return df

    for el in res:
        if "-" in el["alleles"]:
            continue
        if len(el["alleles"]) > 2:
            ref_allele = el["alleles"][0]
            alt_alleles = el["alleles"][1:]
            for allele in alt_alleles:
                new_row = pd.DataFrame(
                    {"ensembl_gene_id": [ens_gene_id],
                     "gene_name": [gene_name],
                     "chromosome": [el["seq_region_name"]],
                     "ref_allele": [ref_allele], "start_pos": [el["start"]],
                     "alt_allele": [allele],
                     "variant": [create_variant_HGVS(el["seq_region_name"],
                                                       el["start"], ref_allele,
                                                       allele)],
                     "dbsnp_id": [el["id"]]})
                df = df.append(new_row, ignore_index=True)
        elif len(el["alleles"]) == 2:
            new_row = pd.DataFrame(
                {"ensembl_gene_id": [ens_gene_id], "gene_name": [gene_name],
                 "chromosome": [el["seq_region_name"]],
                 "ref_allele": [el["alleles"][0]], "start_pos": [el["start"]],
                 "alt_allele": [el["alleles"][1]],
                 "variant": [create_variant_HGVS(el["seq_region_name"],
                                                   el["start"],
                                                   el["alleles"][0],
                                                   el["alleles"][1])],
                 "dbsnp_id": [el["id"]]})
            df = df.append(new_row, ignore_index=True)
        else:
            pass

    df.drop_duplicates(inplace=True)

    return df


async def get_diseases_from_gene(gene: str,
                                 with_vars: bool = False) -> pd.DataFrame:
    """Retrieve diseases associated to a specific gene, using Ensembl.

    :param str gene: query gene Ensembl ID or common name

    :param bool with_vars: True to also retrieve phenotypes associated to
        variants of the gene

    :return: pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
        "disease_name", "phenotype_ids"]) or
        pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "variation",
        "disease_name", "phenotype_ids"]) if with_vars=True
    """
    if gene.startswith("MT-"):
        gene = gene.upper().split("-")[1]

    server = "https://rest.ensembl.org"
    ext = "/phenotype/gene/homo_sapiens/{}?include_associated={}".format(
        gene, int(with_vars)
    )
    res = await get_json_request(server + ext)

    if with_vars:
        df = pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "variation",
                                   "disease_name", "phenotype_ids"])
    else:
        df = pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                   "disease_name", "phenotype_ids"])

    if gene.startswith("ENSG"):  # input gene is actually an Ensembl ID
        try:
            ensembl_gene_id = gene
            gene_name = Mitocarta.query.filter(
                Mitocarta.ensembl_id == ensembl_gene_id
            ).first().gene_symbol
        except AttributeError:  # gene not in Mitocarta
            return df
    else:  # input gene is actually a gene name
        try:
            gene_name = gene
            ensembl_gene_id = Mitocarta.query.filter(
                Mitocarta.gene_symbol == gene_name
            ).first().ensembl_id
        except AttributeError:  # gene not in Mitocarta
            return df

    for el in res:
        if el == "error":
            return df

        if with_vars:
            if "attributes" in el \
                    and "associated_gene" in el["attributes"].keys() \
                    and "Variation" in el and "ontology_accessions" in el:
                row = pd.DataFrame(
                    {"ensembl_gene_id": ensembl_gene_id,
                     "gene_name": [gene_name],
                     "variation": [el["Variation"]],
                     "disease_name": [el["description"]],
                     "phenotype_ids": [el["ontology_accessions"]]}
                )
        else:
            if "ontology_accessions" in el:
                row = pd.DataFrame(
                    {"ensembl_gene_id": [ensembl_gene_id],
                     "gene_name": gene_name,
                     "disease_name": [el["description"]],
                     "phenotype_ids": [el["ontology_accessions"]]}
                )
        try:
            df = df.append(row, ignore_index=True)
        except UnboundLocalError:
            return df
    df.drop_duplicates(subset="disease_name", inplace=True)

    return df


async def json_from_gene(gene_input: str) -> dict:
    """Create the final json structure from gene data.

    :param str gene_input: Ensemble gene ID to use for the queries

    :return: dict json("variants": [variants list],
        "diseases": [diseases list])
    """
    gene_name = Mitocarta.query.filter(
        Mitocarta.ensembl_id == gene_input
    ).first().gene_symbol
    gene_df = pd.DataFrame({"ensembl_gene_id": [gene_input],
                            "gene_name": [gene_name]})
    gene_json = json.loads(gene_df.to_json(orient="records"))
    vars_df = await get_vars_from_gene(gene_input)  # we might want to add phenotypes --> get_diseases_from_dbsnp()
    vars_df["disease_name"] = ""
    vars_df["umls_disease_id"] = ""
    vars_df["disease_id"] = ""

    # add disease name and UMLS ID to variants where available
    var_to_disease = VarDiseaseAss.query.filter(
        VarDiseaseAss.dbsnp_id.in_(vars_df.dbsnp_id)
    ).all()
    for idx in vars_df.index:
        for item in var_to_disease:
            if vars_df.at[idx, "dbsnp_id"] == item.dbsnp_id:
                vars_df.at[idx, "disease_name"] = item.disease_name
                vars_df.at[idx, "umls_disease_id"] = item.umls_disease_id
                vars_df.at[idx, "disease_id"] = umls_to_disease_id(item.umls_disease_id)

    vars_df.drop(["ref_allele", "start_pos", "alt_allele"], axis=1,
                 inplace=True)
    vars_json = json.loads(vars_df.to_json(orient="records"))

    disease_df = await get_diseases_from_gene(gene_input)
    entrez_gene = await ensembl_gene_id_to_entrez(gene_input)
    entrez_gene_id = entrez_gene.entrez_gene_id.values[0]
    gene_to_disease = GeneDiseaseAss.query.filter(
        GeneDiseaseAss.gene_symbol == entrez_gene_id
    ).all()

    if disease_df.shape[0] == 0 and len(gene_to_disease) > 0:
        disease_df = pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                           "disease_id", "disease_name",
                                           "umls_disease_id"])
        for el in gene_to_disease:
            disease_df = disease_df.append(pd.DataFrame(
                {"ensembl_gene_id": [gene_input],
                 "gene_name": [el.gene_symbol],
                 "disease_id": [umls_to_disease_id(el.umls_disease_id)],
                 "disease_name": [el.disease_name],
                 "umls_disease_id": [el.umls_disease_id]}),
                ignore_index=True)
        disease_df["phenotype_ids"] = ""
        disease_df.drop_duplicates(subset="disease_name", inplace=True)
    else:
        disease_df["umls_disease_id"] = ""
        disease_df["disease_id"] = ""
        for idx in disease_df.index:
            for item in gene_to_disease:
                if fuzz.token_sort_ratio(
                        disease_df.at[idx, "disease_name"].capitalize(),
                        item.disease_name.capitalize()
                ) >= 90:
                    disease_df.at[idx, "disease_name"] = item.disease_name
                    disease_df.at[idx, "umls_disease_id"] = item.umls_disease_id
                    disease_df.at[idx, "disease_id"] = umls_to_disease_id(item.umls_disease_id)
                    break

    disease_json = json.loads(disease_df.to_json(orient="records"))
    phenos_df = pd.DataFrame(columns=["phenotype_id", "phenotype_name"])
    for el in disease_json:
        el["phenotype_names"] = []
        for pheno in el["phenotype_ids"]:
            pheno_name = await pheno_id_to_term(pheno)
            el["phenotype_names"].append(pheno_name)
            phenos_df = phenos_df.append(pd.DataFrame(
                {"phenotype_id": [pheno], "phenotype_name": [pheno_name]}
            ), ignore_index=True)
        el["disease_id"] = disease_name_to_id(el["disease_name"])

    phenos_json = json.loads(phenos_df.to_json(orient="records"))

    final_json = dict()
    final_json["variants"] = vars_json
    final_json["genes"] = gene_json
    final_json["diseases"] = disease_json
    final_json["phenotypes"] = phenos_json

    return final_json


def network_from_gene_json(final_json: dict) -> dict:
    """Create the required nodes and edges dictionaries to build the network
    from gene data.

    :param dict final_json: output from json_from_gene()

    :return: dict("nodes": [nodes list], "edges": [edges list])
    """
    var_json = final_json["variants"]
    dis_json = final_json["diseases"]

    v_variants = []
    v_variants.extend([el["variant"] for el in var_json])
    v_variants = set(v_variants)
    v_variants.discard("")

    genes = []
    genes.extend([el["gene_name"] for el in var_json])
    genes = set(genes)
    genes.discard("")

    v_diseases = []
    v_diseases.extend([el["disease_name"]
                       for el in var_json if el["disease_name"] != ""])
    v_diseases = set(v_diseases)
    v_diseases.discard("")

    d_diseases = []
    d_diseases.extend([el["disease_name"] for el in dis_json])
    d_diseases = set(d_diseases)
    d_diseases.discard("")

    d_phenotypes = []
    for el in dis_json:
        d_phenotypes.extend(el["phenotype_names"])
    d_phenotypes = set(d_phenotypes)
    d_phenotypes.discard("")
    d_phenotypes.difference_update(d_diseases)

    nodes = []
    edges = []
    ids = 0
    id_dict = {}
    for el in genes:
        ids += 1
        node = Node("g", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el, "color": {"background": "#739E82",
                                                        "border": "#5F826B"}})
    for el in v_variants:
        ids += 1
        node = Node("v", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el, "color": {"background": "#F9CF45",
                                                        "border": "#CCAA39"}})
    for el in v_diseases:
        ids += 1
        node = Node("d", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                        "border": "#B06A57"}})
    for el in d_diseases:
        ids += 1
        node = Node("d", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                        "border": "#B06A57"}})
    for el in d_phenotypes:
        ids += 1
        node = Node("p", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el, "color": {"background": "#93B5C6",
                                                        "border": "#7995A3"}})

    connected_nodes = set()
    vars_set = set()  # (variant, dbsnp_id, gene_name)
    gene_set = set()  # (gene_name, ensembl_gene_id)
    dise_set = set()  # (disease_name, umls_disease_id)
    phen_set = set()  # (phenotype_name, phenotype_id)

    for el in var_json:
        # gene to variants
        try:
            edges.append({
                "from": id_dict[hash("g" + el["gene_name"])],
                "to": id_dict[hash("v" + el["variant"])]
            })
            connected_nodes.add(id_dict[hash("g" + el["gene_name"])])
            gene_set.add((el["gene_name"], el["ensembl_gene_id"]))
            connected_nodes.add(id_dict[hash("v" + el["variant"])])
            vars_set.add((el["variant"], el["dbsnp_id"], el["gene_name"]))
            # variant to diseases
            if el["disease_name"] != "":
                edges.append({
                    "from": id_dict[hash("v" + el["variant"])],
                    "to": id_dict[hash("d" + el["disease_name"])]
                })
                connected_nodes.add(id_dict[hash("v" + el["variant"])])
                vars_set.add((el["variant"], el["dbsnp_id"], el["gene_name"]))
                connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
                dise_set.add((el["disease_name"], el["disease_id"]))
        except:
            pass
    for el in dis_json:
        # gene to diseases
        try:
            edges.append({
                "from": id_dict[hash("g" + el["gene_name"])],
                "to": id_dict[hash("d" + el["disease_name"])]
            })
            connected_nodes.add(id_dict[hash("g" + el["gene_name"])])
            gene_set.add((el["gene_name"], el["ensembl_gene_id"]))
            connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
            dise_set.add((el["disease_name"], el["disease_id"]))
        except:
            pass
        # disease to phenotypes
        for n, pheno in enumerate(el["phenotype_names"]):
            try:
                if id_dict[hash("d" + el["disease_name"])] != id_dict[hash("p" + pheno)]:
                    edges.append({
                        "from": id_dict[hash("d" + el["disease_name"])],
                        "to": id_dict[hash("p" + pheno)]
                    })
                    connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
                    dise_set.add((el["disease_name"], el["disease_id"]))
                    connected_nodes.add(id_dict[hash("p" + pheno)])
                    phen_set.add((pheno, el["phenotype_ids"][n]))
            except KeyError:
                pass

    # delete orphan nodes
    true_nodes = [el for el in nodes if el["id"] in connected_nodes]

    return {"nodes": true_nodes, "edges": edges, "variants": vars_set,
            "genes": gene_set, "diseases": dise_set, "phenotypes": phen_set}


# FROM PHENOTYPE #


def get_genes_from_phenotype(phenotype: str) -> pd.DataFrame:
    """Retrieve genes related to a phenotype, using Ensembl.

    :param str phenotype: accession id of the phenotype to search for

    :return: pd.DataFrame(columns=["gene_name", "ensembl_gene_id",
        "phenotype_id", "phenotype_name"])
    """
    df = pd.DataFrame(columns=["gene_name", "ensembl_gene_id", "phenotype_id",
                               "phenotype_name"])

    dis_maps = DiseaseMappings.query.filter(
        DiseaseMappings.disease_id == phenotype
    ).first()
    try:
        pheno_umls = dis_maps.umls_disease_id
        pheno_name = dis_maps.disease_name
    except AttributeError:
        return df
    rel_genes = GeneDiseaseAss.query.filter(
        GeneDiseaseAss.umls_disease_id == pheno_umls
    ).all()
    if len(rel_genes) > 0:
        for el in rel_genes:
            ensembl_id = Mitocarta.query.filter(
                Mitocarta.gene_symbol == el.gene_symbol
            ).first().ensembl_id
            row = pd.DataFrame({"gene_name": [el.gene_symbol],
                                "ensembl_gene_id": [ensembl_id],
                                "phenotype_id": [phenotype],
                                "phenotype_name": [pheno_name]})
            df = df.append(row, ignore_index=True)
    else:
        return df

    df.drop_duplicates(inplace=True)

    return df


async def get_vars_from_phenotype(phenotype: str) -> pd.DataFrame:
    """Retrieve variants related to a phenotype, exploiting the
    get_genes_from_phenotype() and get_vars_from_gene_name() functions.

    :param str phenotype: accession id of the phenotype to search for

    :return: pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
        "chromosome", "ref_allele", "start_pos", "alt_allele",
        "phenotype_name", "phenotype_id"])
    """
    dis_maps = DiseaseMappings.query.filter(
        DiseaseMappings.disease_id == phenotype
    ).first()
    if dis_maps is None:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                     "chromosome", "ref_allele", "start_pos",
                                     "alt_allele", "phenotype_name",
                                     "phenotype_id"])
    pheno_umls = dis_maps.umls_disease_id
    pheno_name = dis_maps.disease_name
    vars_maps = VarDiseaseAss.query.filter(
        VarDiseaseAss.umls_disease_id == pheno_umls
    ).all()

    if len(vars_maps) == 0:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                     "chromosome", "ref_allele", "start_pos",
                                     "alt_allele", "phenotype_name",
                                     "phenotype_id"])

    res = await apy.aquery(attributes=["chr_name", "chrom_start",
                                       "consequence_allele_string",
                                       "ensembl_gene_stable_id",
                                       "associated_gene", "refsnp_id",
                                       "phenotype_description"],
                           filters={"snp_filter": [el.dbsnp_id for el in vars_maps]},
                           dataset="hsapiens_snp")
    res.rename({"Chromosome/scaffold name": "chromosome",
                "Variant name": "dbsnp_id",
                "Chromosome/scaffold position start (bp)": "start_pos",
                "Associated gene with phenotype": "gene_name",
                "Consequence specific allele": "ref/alt allele",
                "Gene stable ID": "ensembl_gene_id",
                "Phenotype description": "phenotype_name"}, axis=1,
               inplace=True)
    res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
    res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
    res = res[res["alt_allele"] != "HGMD_MUTATION"]
    # removing weird non-standard chroms
    res = res[res["chromosome"].astype(str).isin(["1", "2", "3", "4", "5", "6",
                                                  "7", "8", "9", "10", "11",
                                                  "12", "13", "14", "15", "16",
                                                  "17", "18", "19", "20", "21",
                                                  "22", "X", "Y", "M", "MT"])]

    variants = [create_variant_HGVS(el.chromosome, el.start_pos,
                                      el.ref_allele, el.alt_allele)
                for el in res.itertuples()]
    res["variant"] = variants
    res["phenotype_name"] = pheno_name
    res["phenotype_id"] = phenotype
    res.drop(["ref/alt allele"], axis=1, inplace=True)
    res["gene_name"] = res["gene_name"].apply(
        lambda x: x.split("-")[1]
            if type(x) == str and x.startswith("MT-") else x)
    res.drop_duplicates(inplace=True)
    res = res[res["variant"] != "chr_:_>_"]
    if res["gene_name"].dtype == str:
        res = res[
            (res["gene_name"] != "intergenic") &
            (res["gene_name"] != "Intergenic") &
            (res["gene_name"].notnull()) &
            (~res["gene_name"].str.contains(",", na=False))
            ]
    else:  # gene_name is not found, so it is a float64 NaN
        genes = []
        for el in res["ensembl_gene_id"]:
            gene_name = ""
            if type(el) == str:
                try:
                    gene_name_entr = await ensembl_gene_id_to_entrez(el)
                    gene_name = gene_name_entr.gene_name.values[0]
                except IndexError:
                    if el.startswith("LRG"):
                        gene_name = el
                genes.append(gene_name)
        res["gene_name"] = genes

    return res


def get_diseases_from_phenotype(phenotype: str) -> pd.DataFrame:
    """Retrieve diseases related to a phenotype, exploiting the
    get_genes_from_phenotype() and get_diseases_from_gene_name() functions.

    :param str phenotype: accession id of the phenotype to search for

    :return: pd.DataFrame(columns=["pheno_id", "pheno_name", "disease_name",
        "disease_id", "umls_disease_id"])
    """
    df = pd.DataFrame(columns=["phenotype_id", "phenotype_name",
                               "disease_name", "disease_id", "umls_disease_id"])
    hpo_dis = HpoDisGenePhen.query.filter(
        HpoDisGenePhen.hpo_id == phenotype
    ).all()
    try:
        pheno_name = DiseaseMappings.query.filter(
            DiseaseMappings.disease_id == phenotype
        ).first().disease_name
    except AttributeError:
        return df

    if hpo_dis:
        for el in hpo_dis:
            dis_id = el.disease_id.split(":")[1]
            dis_name = Omim.query.filter(Omim.mim_number == dis_id).first()
            if dis_name is not None:
                dis_mim_name = dis_name.mim_name
                row = pd.DataFrame(
                    {"phenotype_id": [el.hpo_id],
                     "phenotype_name": [pheno_name],
                     "disease_name": [dis_mim_name],
                     "disease_id": [el.disease_id],
                     "umls_disease_id": [get_umls_from_disease_id(el.disease_id)]}
                )
                df = df.append(row, ignore_index=True)

    df.drop_duplicates(inplace=True)

    return df


async def json_from_phenotype(pheno_input: str) -> dict:
    """Create the final json structure from phenotype data.

    :param str pheno_input: phenotype ID to use for the queries

    :return: dict json("phenotype": phenotype name,
        "variants": [variants list], "genes": [genes list],
        "diseases": [diseases list])
    """
    vars_df = await get_vars_from_phenotype(pheno_input)
    gene_df = get_genes_from_phenotype(pheno_input)
    disease_df = get_diseases_from_phenotype(pheno_input)

    # Variants
    vars_df.drop(["start_pos", "ref_allele", "alt_allele"], axis=1,
                 inplace=True)
    vars_json = json.loads(vars_df.to_json(orient="records"))
    # Genes
    gene_json = json.loads(gene_df.to_json(orient="records"))
    # Diseases
    disease_json = json.loads(disease_df.to_json(orient="records"))

    final_json = dict()
    pheno_name = DiseaseMappings.query.filter(
        DiseaseMappings.disease_id == pheno_input
    ).first().disease_name
    final_json["phenotypes"] = [{"phenotype_name": pheno_name}]
    final_json["variants"] = vars_json
    final_json["genes"] = gene_json
    final_json["diseases"] = disease_json

    return final_json


def network_from_phenotype_json(final_json: dict) -> dict:
    """Create the required nodes and edges dictionaries to build the network
    from phenotype data.

    :param dict final_json: output from json_from_phenotype()

    :return: dict("nodes": [nodes list], "edges": [edges list])
    """
    phenotype = final_json["phenotypes"]
    var_json = final_json["variants"]
    gen_json = final_json["genes"]
    dis_json = final_json["diseases"]

    variants = []
    variants.extend([el["variant"] for el in var_json])
    variants = set(variants)
    variants.discard("")

    v_genes = []
    v_genes.extend([el["gene_name"] for el in var_json])
    v_genes = set(v_genes)
    v_genes.discard("")

    genes = []
    genes.extend([el["gene_name"] for el in gen_json])
    genes = set(genes)
    genes.discard("")

    diseases = []
    diseases.extend([el["disease_name"] for el in dis_json])
    diseases = set(diseases)
    diseases.discard("")

    nodes = []
    edges = []
    ids = 0
    id_dict = {}

    try:
        ids += 1
        node = Node("p", ids, phenotype[0].get("phenotype_name"))
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": phenotype[0].get("phenotype_name"),
                      "color": {"background": "#93B5C6", "border": "#7995A3"}})
    except:
        pass

    for el in variants:
        ids += 1
        node = Node("v", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el,
                      "color": {"background": "#F9CF45", "border": "#CCAA39"}})
    for el in v_genes:
        ids += 1
        node = Node("g", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el,
                      "color": {"background": "#739E82", "border": "#5F826B"}})
    for el in genes:
        ids += 1
        node = Node("g", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el,
                      "color": {"background": "#739E82", "border": "#5F826B"}})
    for el in diseases:
        ids += 1
        node = Node("d", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el,
                      "color": {"background": "#D7816A", "border": "#B06A57"}})

    connected_nodes = set()
    vars_set = set()  # (variant, dbsnp_id, gene_name)
    gene_set = set()  # (gene_name, ensembl_gene_id)
    dise_set = set()  # (disease_name, umls_disease_id)
    phen_set = set()  # (phenotype_name, phenotype_id)

    for el in var_json:
        # phenotype to variants
        try:
            edges.append({
                "from": id_dict[hash("p" + el["phenotype_name"])],
                "to": id_dict[hash("v" + el["variant"])]
            })
            connected_nodes.add(id_dict[hash("p" + el["phenotype_name"])])
            phen_set.add((el["phenotype_name"], el["phenotype_id"]))
            connected_nodes.add(id_dict[hash("v" + el["variant"])])
            vars_set.add((el["variant"], el["dbsnp_id"], el["gene_name"]))
            # variant to gene
            edges.append({
                "from": id_dict[hash("v" + el["variant"])],
                "to": id_dict[hash("g" + el["gene_name"])]
            })
            connected_nodes.add(id_dict[hash("v" + el["variant"])])
            vars_set.add((el["variant"], el["dbsnp_id"], el["gene_name"]))
            connected_nodes.add(id_dict[hash("g" + el["gene_name"])])
            gene_set.add((el["gene_name"], el["ensembl_gene_id"]))
        except:
            pass
    for el in gen_json:
        # phenotype to genes
        try:
            edges.append({
                "from": id_dict[hash("p" + el["phenotype_name"])],
                "to": id_dict[hash("g" + el["gene_name"])]
            })
            connected_nodes.add(id_dict[hash("p" + el["phenotype_name"])])
            phen_set.add((el["phenotype_name"], el["phenotype_id"]))
            connected_nodes.add(id_dict[hash("g" + el["gene_name"])])
            gene_set.add((el["gene_name"], el["ensembl_gene_id"]))
        except:
            pass
    for el in dis_json:
        # phenotype to diseases
        try:
            edges.append({
                "from": id_dict[hash("p" + el["phenotype_name"])],
                "to": id_dict[hash("d" + el["disease_name"])]
            })
            connected_nodes.add(id_dict[hash("p" + el["phenotype_name"])])
            phen_set.add((el["phenotype_name"], el["phenotype_id"]))
            connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
            dise_set.add((el["disease_name"], el["disease_id"]))
        except:
            pass

    # delete orphan nodes
    true_nodes = [el for el in nodes if el["id"] in connected_nodes]

    return {"nodes": true_nodes, "edges": edges, "variants": vars_set,
            "genes": gene_set, "diseases": dise_set, "phenotypes": phen_set}


# FROM DISEASE #


def get_umls_from_disease_id(disease_id: str) -> str:
    """Convert the general disease ID to the standard UMLS ID.

    :param str disease_id: disease ID starting with "DO", "MSH", "NCI",
        "OMIM", "ORDO" or "ICD9CM"

    :return: str the correspondent UMLS ID
    """
    dis_vocab, dis_num = disease_id.split(":")
    dis_onts = ["DO", "MSH", "NCI", "OMIM", "ORDO", "ICD9CM"]
    if dis_vocab in dis_onts:
        q = DiseaseMappings.query.filter(
            DiseaseMappings.vocabulary == dis_vocab,
            DiseaseMappings.disease_id == dis_num
        ).first()
    elif dis_vocab == "HPO":
        q = DiseaseMappings.query.filter(
            DiseaseMappings.disease_id == disease_id
        ).first()
    elif dis_vocab == "EFO":
        q = DiseaseMappings.query.filter(
            DiseaseMappings.disease_id == "EFO_{}".format(dis_num)
        ).first()
    else:
        return ""

    if q:
        disease_umls = q.umls_disease_id
        return disease_umls
    return ""


def get_genes_from_disease_id(disease_id: str) -> pd.DataFrame:
    """Retrieve genes involved in a specific disease, using the GeneDiseaseAss
    table from the db.

    :param str disease_id: disease ID to look for

    :return: pd.DataFrame(columns=["umls_disease_id", "disease_name",
        "disease_id", "entrez_gene_id", "gene_name", "ensembl_gene_id",
        "ass_score"])
    """
    dis_umls = get_umls_from_disease_id(disease_id)
    ass_genes = GeneDiseaseAss.query.filter(
        GeneDiseaseAss.umls_disease_id == dis_umls
    ).all()
    df = pd.DataFrame(columns=["umls_disease_id", "disease_name",
                               "disease_id", "entrez_gene_id", "gene_name",
                               "ensembl_gene_id", "ass_score"])
    if len(ass_genes) > 0:
        for el in ass_genes:
            ens_gene_id = ""
            mitoq = Mitocarta.query.filter(
                Mitocarta.gene_symbol == el.gene_symbol
            ).first()
            if mitoq is not None:
                ens_gene_id = mitoq.ensembl_id
            row = pd.DataFrame({"disease_id": [disease_id],
                                "umls_disease_id": [dis_umls],
                                "disease_name": [el.disease_name],
                                "entrez_gene_id": [el.entrez_gene_id],
                                "gene_name": [el.gene_symbol],
                                "ensembl_gene_id": [ens_gene_id],
                                "ass_score": [el.score]})
            df = df.append(row, ignore_index=True)
    else:
        ass_genes = HpoDisGenePhen.query.filter(
            HpoDisGenePhen.disease_id == disease_id
        ).all()
        for el in ass_genes:
            ens_gene_id = ""
            mitoq = Mitocarta.query.filter(
                Mitocarta.gene_symbol == el.gene_symbol
            ).first()
            if mitoq is not None:
                ens_gene_id = mitoq.ensembl_id
            row = pd.DataFrame({"disease_id": [disease_id],
                                "umls_disease_id": [dis_umls],
                                "disease_name": [disease_id_to_name(disease_id)],
                                "entrez_gene_id": [el.entrez_gene_id],
                                "gene_name": [el.gene_symbol],
                                "ensembl_gene_id": [ens_gene_id],
                                "ass_score": [0.0]})
            df = df.append(row, ignore_index=True)

    df.drop_duplicates(inplace=True)

    return df


async def get_vars_from_disease_id(disease_id: str) -> pd.DataFrame:
    """Retrieve variants associated to a specific disease, getting their
    dbSNP ID and then the actual variant string exploiting Ensembl API.

    :param str disease_id: disease ID to look for

    :return: pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "dbsnp_id",
        "variant", "umls_disease_id", "disease_name", "disease_id"])
    """
    dis_umls = get_umls_from_disease_id(disease_id)
    ass_vars = VarDiseaseAss.query.filter(
        VarDiseaseAss.umls_disease_id == dis_umls
    ).all()
    dis_maps = DiseaseMappings.query.filter(
        DiseaseMappings.umls_disease_id == dis_umls
    ).first()
    dis_name = ""
    if dis_maps is not None:
        dis_name = dis_maps.disease_name
    if len(ass_vars) == 0:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                     "dbsnp_id", "variant", "umls_disease_id",
                                     "disease_name", "disease_id"])

    res = await apy.aquery(attributes=["chr_name", "chrom_start",
                                       "consequence_allele_string",
                                       "ensembl_gene_stable_id",
                                       "associated_gene", "refsnp_id"],
                           filters={"snp_filter": [el.dbsnp_id for el in ass_vars]},
                           dataset="hsapiens_snp")
    res.rename({"Chromosome/scaffold name": "chromosome",
                "Variant name": "dbsnp_id",
                "Chromosome/scaffold position start (bp)": "start_pos",
                "Associated gene with phenotype": "gene_name",
                "Consequence specific allele": "ref/alt allele",
                "Gene stable ID": "ensembl_gene_id"}, axis=1, inplace=True)
    if res.shape[0] == 0:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name",
                                     "dbsnp_id", "variant", "umls_disease_id",
                                     "disease_name", "disease_id"])
    res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
    res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
    res = res[res["alt_allele"] != "HGMD_MUTATION"]
    res["variant"] = [create_variant_HGVS(el.chromosome, el.start_pos,
                                            el.ref_allele, el.alt_allele)
                      for el in res.itertuples()]
    res["umls_disease_id"] = dis_umls
    res["disease_name"] = dis_name
    res.drop(["ref/alt allele", "ref_allele", "alt_allele", "chromosome",
              "start_pos"], axis=1, inplace=True)
    res["gene_name"] = res["gene_name"].apply(
        lambda x: x.split("-")[1]
            if type(x) == str and x.startswith("MT-") else x)
    gene_dict = dict(zip(res["ensembl_gene_id"], res["gene_name"]))
    res["gene_name"] = res.apply(
        lambda row: gene_dict[row["ensembl_gene_id"]]
        if type(row["gene_name"]) != str else row["gene_name"],
        axis=1
    )
    if res["gene_name"].dtype == str:
        res = res[(res["gene_name"] != "intergenic") &
                  (res["gene_name"] != "Intergenic") &
                  (res["gene_name"].notnull()) &
                  (~res["gene_name"].str.contains(",", na=False))]
    res["disease_id"] = disease_id
    res.drop_duplicates(subset="variant", inplace=True)

    return res


def get_phenos_from_disease_id(disease_id: str) -> pd.DataFrame:
    """Retrieve phenotypes associated to a specific disease, using the
    HpoDisGenePhen table from the db.

    :param str disease_id: disease ID to look for

    :return: pd.DataFrame(columns=["umls_disease_id", "disease_name",
        "disease_id", "phenotype_id", "phenotype_name"]
    """
    df = pd.DataFrame(columns=["umls_disease_id", "disease_name", "disease_id",
                               "phenotype_id", "phenotype_name"])
    ass_phenos = HpoDisGenePhen.query.filter(
        HpoDisGenePhen.disease_id == disease_id
    ).all()
    try:
        disease_name = Diseases.query.filter(
            Diseases.disease_id == disease_id
        ).first().disease_name
    except AttributeError:
        return df
    dis_maps = DiseaseMappings.query.filter(
        DiseaseMappings.umls_disease_id == get_umls_from_disease_id(disease_id)
    ).first()
    disease_umls = ""
    if dis_maps is not None:
        disease_name = dis_maps.disease_name
        disease_umls = dis_maps.umls_disease_id

    if ass_phenos:
        for el in ass_phenos:
            row = pd.DataFrame({"umls_disease_id": disease_umls,
                                "disease_name": [disease_name],
                                "disease_id": [disease_id],
                                "phenotype_id": [el.hpo_id],
                                "phenotype_name": [el.hpo_term_name]})
            df = df.append(row, ignore_index=True)

    df.drop_duplicates(inplace=True)

    return df


async def json_from_disease(disease_input: str) -> dict:
    """Create the final json structure from disease data.

    :param str disease_input: disease ID to use for the queries

    :return: dict json("diseases": disease name,
        "phenotype": [phenotypes list], "variants": [variants list],
        "genes": [genes list])
    """
    vars_df = await get_vars_from_disease_id(disease_input)
    gene_df = get_genes_from_disease_id(disease_input)
    pheno_df = get_phenos_from_disease_id(disease_input)

    # Variants
    vars_json = json.loads(vars_df.to_json(orient="records"))
    # Genes
    gene_df.drop(["entrez_gene_id", "ass_score"], axis=1, inplace=True)  # do not need these for now
    gene_json = json.loads(gene_df.to_json(orient="records"))
    # Phenotypes
    pheno_json = json.loads(pheno_df.to_json(orient="records"))

    final_json = dict()
    disease_name = Diseases.query.filter(
        Diseases.disease_id == disease_input
    ).first().disease_name
    dis_maps = DiseaseMappings.query.filter(
        DiseaseMappings.umls_disease_id == get_umls_from_disease_id(disease_input)
    ).first()
    if dis_maps is not None:
        disease_name = dis_maps.disease_name

    final_json["diseases"] = [{"disease_name": disease_name}]
    final_json["phenotypes"] = pheno_json
    final_json["variants"] = vars_json
    final_json["genes"] = gene_json

    return final_json


def network_from_disease_json(final_json: dict) -> dict:
    """Create the required nodes and edges dictionaries to build the network
    from disease data.

    :param dict final_json: output from json_from_disease()

    :return: dict("nodes": [nodes list], "edges": [edges list])
    """
    phen_json = final_json["phenotypes"]
    var_json = final_json["variants"]
    gen_json = final_json["genes"]
    disease = final_json["diseases"]
    # dis_json = final_json["diseases"]

    variants = []
    variants.extend([el["variant"] for el in var_json])
    variants = set(variants)
    variants.discard("")

    v_genes = []
    v_genes.extend([el["gene_name"] for el in var_json])
    v_genes = set(v_genes)
    v_genes.discard("")

    genes = []
    genes.extend([el["gene_name"] for el in gen_json])
    genes = set(genes)
    genes.discard("")

    phenos = []
    phenos.extend([el["phenotype_name"] for el in phen_json])
    phenos = set(phenos)
    phenos.discard("")

    nodes = []
    edges = []
    ids = 0
    id_dict = {}

    try:
        ids += 1
        node = Node("d", ids, disease[0].get("disease_name"))
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": disease[0].get("disease_name"),
                      "color": {"background": "#D7816A", "border": "#B06A57"}})
    except:
        pass

    for el in variants:
        ids += 1
        node = Node("v", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el,
                      "color": {"background": "#F9CF45", "border": "#CCAA39"}})
    for el in v_genes:
        ids += 1
        if el is not None:
            node = Node("g", ids, el)
            id_dict[node.uid] = ids
            nodes.append({"id": ids, "label": el,
                          "color": {"background": "#739E82", "border": "#5F826B"}})
    for el in genes:
        ids += 1
        node = Node("g", ids, el)
        id_dict[node.uid] = ids
        nodes.append({"id": ids, "label": el,
                      "color": {"background": "#739E82", "border": "#5F826B"}})
    for el in phenos:
        if el != disease[0].get("disease_name"):
            ids += 1
            node = Node("p", ids, el)
            id_dict[node.uid] = ids
            nodes.append({"id": ids, "label": el,
                          "color": {"background": "#93B5C6",
                                    "border": "#7995A3"}})

    connected_nodes = set()
    vars_set = set()  # (variant, dbsnp_id, gene_name)
    gene_set = set()  # (gene_name, ensembl_gene_id)
    dise_set = set()  # (disease_name, umls_disease_id)
    phen_set = set()  # (phenotype_name, phenotype_id)

    for el in var_json:
        # disease to variants
        try:
            edges.append({
                "from": id_dict[hash("d" + el["disease_name"])],
                "to": id_dict[hash("v" + el["variant"])]
            })
            connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
            dise_set.add((el["disease_name"], el["disease_id"]))
            connected_nodes.add(id_dict[hash("v" + el["variant"])])
            vars_set.add((el["variant"], el["dbsnp_id"], el["gene_name"]))
        except:
            pass
        # variants to genes
        try:
            edges.append({
                "from": id_dict[hash("v" + el["variant"])],
                "to": id_dict[hash("g" + el["gene_name"])]
            })
            connected_nodes.add(id_dict[hash("v" + el["variant"])])
            vars_set.add((el["variant"], el["dbsnp_id"], el["gene_name"]))
            connected_nodes.add(id_dict[hash("g" + el["gene_name"])])
            gene_set.add((el["gene_name"], el["ensembl_gene_id"]))
        except TypeError:
            pass
    for el in gen_json:
        # disease to genes
        try:
            edges.append({
                "from": id_dict[hash("d" + el["disease_name"])],
                "to": id_dict[hash("g" + el["gene_name"])]
            })
            connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
            dise_set.add((el["disease_name"], el["disease_id"]))
            connected_nodes.add(id_dict[hash("g" + el["gene_name"])])
            gene_set.add((el["gene_name"], el["ensembl_gene_id"]))
        except:
            pass
    for el in phen_json:
        # disease to phenotypes
        try:
            if id_dict[hash("p" + el["phenotype_name"])] != "" \
                    and id_dict[hash("d" + el["disease_name"])] != "":
                if id_dict[hash("p" + el["phenotype_name"])] != id_dict[hash("d" + el["disease_name"])]:
                    edges.append({
                        "from": id_dict[hash("d" + el["disease_name"])],
                        "to": id_dict[hash("p" + el["phenotype_name"])]
                    })
                    connected_nodes.add(id_dict[hash("d" + el["disease_name"])])
                    dise_set.add((el["disease_name"], el["disease_id"]))
                    connected_nodes.add(id_dict[hash("p" + el["phenotype_name"])])
                    phen_set.add((el["phenotype_name"], el["phenotype_id"]))
        except:
            pass

    # delete orphan nodes
    true_nodes = [el for el in nodes if el["id"] in connected_nodes]

    return {"nodes": true_nodes, "edges": edges, "variants": vars_set,
            "genes": gene_set, "diseases": dise_set, "phenotypes": phen_set}


def create_dataframes(json_data: dict) -> dict:
    """Create 4 different dataframes for variants, genes, diseases and
    phenotypes from results of json_from_*() functions.

    :param json_data: output of json_from_*() function
    :return: dict of pd.DataFrames
    """
    variants = pd.DataFrame.from_records(json_data["variants"])
    genes = pd.DataFrame.from_records(json_data["genes"])
    diseases = pd.DataFrame.from_records(json_data["diseases"])
    phenotypes = pd.DataFrame.from_records(json_data["phenotypes"])

    return {"variants": variants, "genes": genes,
            "diseases": diseases, "phenotypes": phenotypes}
