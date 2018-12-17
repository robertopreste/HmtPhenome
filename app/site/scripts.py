#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from app.site.models import Mitocarta, Phenotypes, Diseases, Omim, Orphanet, GeneDiseaseAss, VarDiseaseAss, DiseaseMappings
import pandas as pd
import numpy as np
import requests
# import sys
from pybiomart import Server
from fuzzywuzzy import fuzz
import json


def get_genes():
    """
    Retrieve Mitocarta genes from the database and create a dictionary that will be used to
    populate the related dropdown menu in query page.
    :return: dictionary with chromosome as keys and genes as values
    """
    chr_dict = {}
    chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
            "chr20", "chr21", "chr22", "chrX", "chrM"]
    for chrom in chrs:
        lista = set()
        q = Mitocarta.query.filter(Mitocarta.hg_chr == chrom).all()
        for el in q:
            lista.add(el.gene_symbol)
        chr_dict[chrom] = sorted(list(lista))

    return chr_dict


def get_phenos():
    """
    Retrieve phenotypes from the database and create a list that will be used to populate the
    related dropdown menu in the query page.
    :return: list of phenotypes
    """
    pheno_list = []
    q = Phenotypes.query.all()
    for el in q:
        pheno_list.append([el.hpo_term_name, el.hpo_id])

    return pheno_list


def get_diseases():
    """
    Retrieve diseases from the database and create a list that will be used to populate the
    related dropdown menu in the query page.
    :return: list of diseases
    """
    disease_list = []
    q = Diseases.query.all()
    for el in q:
        disease_list.append([el.disease_name, el.disease_id])

    return disease_list


def populate_phenos():
    """
    Add phenotypes data for the autocomplete function in the script.js file, during the update of
    the database.
    :return: JavaScript autocomplete function for phenotypes
    """
    phenos = get_phenos()
    base_string = """
var pheno_compl = document.getElementById("pheno_input"); 
new Awesomplete(pheno_compl, {list: %s});     
    """ % repr(phenos)

    return base_string


def populate_diseases():
    """
    Add diseases data for the autocomplete function in the script.js file, during the update of the
    database.
    :return: JavaScript autocomplete function for diseases
    """
    diseases = get_diseases()
    base_string = """
var disease_compl = document.getElementById("disease_input"); 
new Awesomplete(disease_compl, {list: %s});     
    """ % repr(diseases)

    return base_string


def populate_genes():
    """
    Add genes for each chromosome in the script.js file, during the update of the database, in
    order to populate the related dropdown menu in the query page.
    :return: JavaScript function for genes/chromosomes
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
        "%s|%s", """ % (el, el)
        else:
            base_string += """]; 
    } else if (s1.value == "%s") {
        optionArray = ["A|--All Genes--", """ % chrom

            for el in chr_dict[chrom]:
                base_string += """
        "%s|%s", """ % (el, el)

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


def pheno_name_to_id(pheno_name):
    """
    Retrieve the related HP id from a given phenotype name.
    :param pheno_name: [str] phenotype name
    :return: [list] with the related id(s)
    """
    r = requests.get("https://hpo.jax.org/api/hpo/search?q={}".format(pheno_name),
                     headers={"Content-Type": "application/json"})
    res = r.json()
    ids = []

    if res["termsTotalCount"] == 0:
        pass
    elif res["termsTotalCount"] == 1:
        ids.append(res["terms"][0]["id"])
    else:
        for el in res["terms"]:
            ids.append(el["id"])

    return ids


def pheno_id_to_term(pheno_id):
    """
    Retrieve the common phenotype name from a given ID. It can retrieve the data either from the
    local database (for HP: IDs) or from the web (for EFO: IDs).
    :param pheno_id: [str] HP: or EFO: ID
    :return: [str] with the related common name
    """
    if pheno_id.startswith("HP:"):
        try:
            pheno_name = Phenotypes.query.filter(Phenotypes.hpo_id == pheno_id).first().hpo_term_name
        except AttributeError:
            pheno_name = ""
    elif pheno_id.startswith("Orphanet:"):
        orpha_id = pheno_id.strip("Orphanet:")
        pheno_name = Orphanet.query.filter(Orphanet.orpha_num == orpha_id).first().orpha_name
    elif pheno_id.startswith("EFO:"):
        efo_id = pheno_id.strip("EFO:")
        base_url = "https://www.ebi.ac.uk/ols/api/ontologies/efo/terms?"
        iri_url = "iri=http://www.ebi.ac.uk/efo/EFO_{}".format(efo_id)
        r = requests.get(base_url + iri_url, headers={"Content-Type": "application/json"})
        res = r.json()
        pheno_name = res["_embedded"]["terms"][0]["label"]

    return pheno_name


def disease_id_to_name(disease_id):
    """
    Convert a given disease ID into its common name.
    :param disease_id: [str] query disease ID
    :return: [str] related disease name
    """
    try:
        q = Diseases.query.filter(Diseases.disease_id == disease_id).first()
        return q.disease_name
    except AttributeError:
        try:
            if disease_id.startswith("OMIM:"):
                q = Omim.query.filter(Omim.mim_number == int(disease_id.strip("OMIM:"))).first()
                return q.mim_name
            elif disease_id.startswith("ORPHA:"):
                q = Orphanet.query.filter(Orphanet.orpha_num == int(disease_id.strip("ORPHA:"))).first()
                return q.orpha_name
        except AttributeError:
            return ""


def disease_name_to_id(disease_name):
    """
    Convert a given disease name into the related Omim or Orphanet ID.
    :param disease_name: [str] query disease name
    :return: [str] related disease ID from OMIM or ORPHANET
    """
    try:
        q = Diseases.query.filter(Diseases.disease_name == disease_name).first()
        return q.disease_id
    except AttributeError:
        try:
            q = Omim.query.filter(Omim.mim_name == disease_name).first()
            return "OMIM:{}".format(q.mim_number)  # TODO: check also mim symbol
        except AttributeError:
            try:
                q = Orphanet.query.filter(Orphanet.orpha_name == disease_name).first()
                return "ORPHA:{}".format(q.orpha_num)
            except AttributeError:
                return ""


# FROM VARIANT POSITION #


def get_gene_from_variant(chrom, var_start, var_end=None):
    """
    Retrieve the gene to which the provided variant belongs, using Biomart.
    :param chrom: [str] chromosome name (chr + 1:22, X, Y, M)
    :param var_start: [str, int] variant starting position
    :param var_end: [str, int] variant ending position
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name"]
    """

    chrom = chrom.lstrip("chr").upper()
    if chrom == "M":
        chrom = "MT"
    if var_end is None:
        var_end = var_start

    server = Server(host="http://www.ensembl.org")
    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]

    res = dataset.query(attributes=["ensembl_gene_id", "external_gene_name"],
                        filters={"chromosome_name": chrom, "start": var_start, "end": var_end})

    res.rename({"Gene stable ID": "ensembl_gene_id", "Gene name": "gene_name"},
               axis=1, inplace=True)

    return res


def get_pheno_from_variant(chrom, var_start, var_end=None):
    """
    Retrieve phenotypes associated with a specific variant, using Biomart.
    :param chrom: [str] chromosome name (chr + 1:22, X, Y, M)
    :param var_start: [str, int] variant starting position
    :param var_end: [str, int] variant ending position
    :return: pd.DataFrame with columns ["chromosome", "ref_allele", "start_pos", "alt_allele",
    "phenotype"]
    """

    chrom = chrom.lstrip("chr").upper()
    if chrom == "M":
        chrom = "MT"
    if var_end is None:
        var_end = var_start

    server = Server(host="http://www.ensembl.org")
    dataset = server.marts["ENSEMBL_MART_SNP"].datasets["hsapiens_snp"]

    res = dataset.query(attributes=["chr_name", "chrom_start", "consequence_allele_string",
                                    "phenotype_description"],
                        filters={"chr_name": chrom, "start": var_start, "end": var_end})

    res.rename({"Chromosome/scaffold name": "chromosome",
                "Chromosome/scaffold position start (bp)": "start_pos",
                "Consequence specific allele": "ref/alt allele",
                "Phenotype description": "phenotype"}, axis=1, inplace=True)
    if len(res) > 0:
        res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
        res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
        res = res[res["alt_allele"] != "HGMD_MUTATION"]
        res = res[["chromosome", "ref_allele", "start_pos", "alt_allele", "phenotype"]]
        res.drop_duplicates(subset="phenotype", inplace=True)
    else:
        res = pd.DataFrame(columns=["chromosome", "ref_allele", "start_pos", "alt_allele", "phenotype"])
    res["phenotype_id"] = res["phenotype"].apply(pheno_name_to_id)

    return res


def get_diseases_from_variant(chrom, var_start, var_end=None):
    """
    Retrieve diseases associated with the given variant, exploiting the get_gene_from_variant and
    get_diseases_from_gene_name functions.
    :param chrom: [str] chromosome name (chr + 1:22, X, Y, M)
    :param var_start: [str, int] variant starting position
    :param var_end: [str, int] variant ending position
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "location", "variation",
    "disease", "phenotypes"]
    """
    chrom = chrom.lstrip("chr").upper()
    if chrom == "M":
        chrom = "MT"
    # if var_end is None:
    #     var_end = var_start

    gene_name = get_gene_from_variant(chrom, var_start, var_end)["gene_name"][0]
    diseases = get_diseases_from_gene_name(gene_name, True)

    try:
        if var_end is not None:
            diseases = diseases[(diseases["location"].str.startswith(chrom)) &
                                (diseases["location"].str.contains(str(var_start))) &
                                (diseases["location"].str.contains(var_end))]
        else:
            diseases = diseases[(diseases["location"].str.startswith(chrom)) &
                                (diseases["location"].str.contains(str(var_start)))]
    except KeyError:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "location", "variation",
                                     "disease", "phenotypes"])

    return diseases


def json_from_variant(variant_chr, variant_start, variant_end):
    """
    Create the final json structure from variant data.
    :param disease_df: result of get_diseases_from_variant()
    :param pheno_df: result of get_pheno_from_variant()
    :return:
    """
    disease_df = get_diseases_from_variant(variant_chr, variant_start, variant_end)
    pheno_df = get_pheno_from_variant(variant_chr, variant_start, variant_end)

    df = (disease_df.set_index("disease")
          .join(pheno_df.set_index("phenotype"))
          .reset_index())
    df["variant"] = df["ref_allele"] + df["start_pos"].astype(str) + df["alt_allele"]

    disgen_disease = VarDiseaseAss.query.filter(VarDiseaseAss.dbsnp_id.in_(df.variation)).all()
    for idx in df.index:
        for item in disgen_disease:
            if fuzz.token_sort_ratio(df.at[idx, "disease"].capitalize(),
                                     item.disease_name.capitalize()) >= 90:
                df.at[idx, "disease"] = item.disease_name

    df.drop(["location", "ref_allele", "start_pos", "alt_allele", "phenotype_id"], axis=1,
            inplace=True)
    df.rename({"disease": "disease_name", "variation": "dbsnp_id", "phenotypes": "phenotype_ids"},
              axis=1, inplace=True)

    df_json = json.loads(df.to_json(orient="records"))
    for el in df_json:
        el["phenotype_names"] = []
        for pheno in el["phenotype_ids"]:
            el["phenotype_names"].append(pheno_id_to_term(pheno))
        el["disease_id"] = disease_name_to_id(el["disease_name"])  # TODO: correct this to grep disease id from the db table

    final_json = {}
    final_json["variants"] = df_json

    return final_json


def final_from_variant(gene_df, pheno_df, disease_df):
    """
    Create the final dataframe for variant data.
    :param gene_df: result of get_gene_from_variant()
    :param pheno_df: result of get_pheno_from_variant()
    :param disease_df: result of get_diseases_from_variant()
    :return: pd.DataFrame with columns ["variant", "chromosome", "ensembl_gene_id", "gene_name",
    "variation", "disease_id", "disease", "phenotype_id", "phenotype_name"]
    """
    df = (disease_df.set_index("disease")
          .join(pheno_df.set_index("phenotype"))
          .reset_index())
    df["variant"] = df["ref_allele"] + df["start_pos"].astype(str) + df["alt_allele"]

    final_df = pd.DataFrame(columns=["variant", "chromosome", "ensembl_gene_id", "gene_name",
                                     "variation", "disease_id", "disease", "phenotype_id",
                                     "phenotype_name"])

    for row in df.itertuples():
        disease_id = disease_name_to_id(row.disease)
        for pheno in row.phenotypes:
            new_row = pd.DataFrame({"variant": row.variant, "chromosome": row.chromosome,
                                    "ensembl_gene_id": row.ensembl_gene_id,
                                    "gene_name": row.gene_name, "variation": row.variation,
                                    "disease_id": disease_id, "disease": row.disease,
                                    "phenotype_id": pheno,
                                    "phenotype_name": pheno_id_to_term(pheno)}, index=[row.variant])
            final_df = final_df.append(new_row, ignore_index=True)

    return final_df


def network_from_variant(final_df):
    """
    Create nodes and edges lists for the creation of the network starting from variants.
    :param final_df: final dataframe returned by final_from_variant()
    :return: dictionary with a nodes list and an edges list
    """
    nodes = []
    edges = []
    ids = 0
    id_dict = {}
    variants = final_df["variant"].unique()
    for el in variants:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#F9CF45",
                                                        "border": "#CCAA39"}})
        id_dict[el] = ids
    genes = final_df["gene_name"].unique()
    for el in genes:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#739E82",
                                                        "border": "#5F826B"}})
        id_dict[el] = ids
    diseases = final_df["disease"].unique()
    for el in diseases:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                        "border": "#B06A57"}})
        id_dict[el] = ids
    phenotypes = final_df["phenotype_name"].unique()
    for el in phenotypes:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#93B5C6",
                                                        "border": "#7995A3"}})
        id_dict[el] = ids

    for var in variants:
        gene_names = final_df[final_df["variant"] == var].gene_name.unique()
        dis_names = final_df[final_df["variant"] == var].disease.unique()
        for gene in gene_names:
            edges.append({"from": id_dict[var], "to": id_dict[gene]})
        for dis in dis_names:
            edges.append({"from": id_dict[var], "to": id_dict[dis]})
    for dis in diseases:
        pheno_names = final_df[final_df["disease"] == dis].phenotype_name.unique()
        for pheno in pheno_names:
            edges.append({"from": id_dict[dis], "to": id_dict[pheno]})

    return {"nodes": nodes, "edges": edges}


def network_from_variant_json(final_json):
    """
    Create the required nodes and edges dictionaries to build the network from variant data.
    :param final_json: output from json_from_variant()
    :return: dictionary with nodes and edges for network construction
    """
    var_json = final_json["variants"]
    variants = []
    variants.extend([el["variant"] for el in var_json])
    variants = set(variants)

    genes = []
    genes.extend([el["gene_name"] for el in var_json])
    genes = set(genes)

    diseases = []
    diseases.extend([el["disease_name"] for el in var_json])
    diseases = set(diseases)

    phenotypes = []
    for el in var_json:
        phenotypes.extend(el["phenotype_names"])
    phenotypes = set(phenotypes)

    nodes = []
    edges = []
    ids = 0
    id_dict = {}
    for el in variants:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#F9CF45",
                                                        "border": "#CCAA39"}})
        id_dict[el] = ids
    for el in genes:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#739E82",
                                                        "border": "#5F826B"}})
        id_dict[el] = ids
    for el in diseases:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                        "border": "#B06A57"}})
        id_dict[el] = ids
    for el in phenotypes:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#93B5C6",
                                                        "border": "#7995A3"}})
        id_dict[el] = ids

    # TODO: remove duplicated entries from the edges list

    for el in var_json:
        # variant to gene
        edges.append({"from": id_dict[el["variant"]], "to": id_dict[el["gene_name"]]})
        # variant to diseases
        edges.append({"from": id_dict[el["variant"]], "to": id_dict[el["disease_name"]]})
        for pheno in el["phenotype_names"]:
            # disease to phenotypes
            edges.append({"from": id_dict[el["disease_name"]], "to": id_dict[pheno]})

    return {"nodes": nodes, "edges": edges}


# FROM GENE #


def get_vars_from_gene_name(gene_name):
    """
    Retrieve all variants associated with a specific gene, using Biomart.
    :param gene_name: [str] name of the query gene
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "chromosome", "ref_allele",
    "start_pos", "alt_allele", "variant", "phenotype"]
    """
    # gene_name = gene_name.upper().lstrip("MT-")
    if gene_name.startswith("MT-"):
        gene_name = gene_name.upper().split("-")[1]

    try:
        ens_gene_id = Mitocarta.query.filter(Mitocarta.gene_symbol == gene_name).first().ensembl_id
    except AttributeError:  # gene not in Mitocarta
        return pd.DataFrame()

    server = Server(host="http://www.ensembl.org")
    dataset = server.marts["ENSEMBL_MART_SNP"].datasets["hsapiens_snp"]

    res = dataset.query(attributes=["chr_name", "chrom_start", "consequence_allele_string",
                                    "phenotype_description", "refsnp_id"],
                        filters={"ensembl_gene": ens_gene_id})

    res.rename({"Chromosome/scaffold name": "chromosome", "Variant name": "dbsnp_id",
                "Chromosome/scaffold position start (bp)": "start_pos",
                "Consequence specific allele": "ref/alt allele",
                "Phenotype description": "phenotype"}, axis=1, inplace=True)
    res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
    res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
    res = res[res["alt_allele"] != "HGMD_MUTATION"]
    res["gene_name"] = gene_name
    res["ensembl_gene_id"] = ens_gene_id
    res["variant"] = res["ref_allele"] + res["start_pos"].astype(str) + res["alt_allele"]

    # TODO: add also variants found through the VarDiseaseAss table (then retrieve variant from dbsnp id)

    res = res[["ensembl_gene_id", "gene_name", "chromosome", "ref_allele", "start_pos", "alt_allele",
               "variant", "dbsnp_id", "phenotype"]]
    res = res[res["phenotype"].notnull()]

    return res


def get_vars_from_gene_id(ens_gene_id):
    """
    Retrieve all variants associated with a specific Ensembl gene ID, using Biomart.
    :param ens_gene_id: [str] query Ensembl gene ID
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "chromosome", "ref_allele",
    "start_pos", "alt_allele", "phenotype"]
    """
    try:
        gene_name = Mitocarta.query.filter(Mitocarta.ensembl_id == ens_gene_id).first().gene_symbol
    except AttributeError:  # gene not in Mitocarta
        return pd.DataFrame()

    server = Server(host="http://www.ensembl.org")
    dataset = server.marts["ENSEMBL_MART_SNP"].datasets["hsapiens_snp"]

    res = dataset.query(attributes=["chr_name", "chrom_start", "consequence_allele_string",
                                    "phenotype_description"],
                        filters={"ensembl_gene": ens_gene_id})

    res.rename({"Chromosome/scaffold name": "chromosome",
                "Chromosome/scaffold position start (bp)": "start_pos",
                "Consequence specific allele": "ref/alt allele",
                "Phenotype description": "phenotype"}, axis=1, inplace=True)
    res["ref_allele"] = res["ref/alt allele"].str.split("/", expand=True)[0]
    res["alt_allele"] = res["ref/alt allele"].str.split("/", expand=True)[1]
    res = res[res["alt_allele"] != "HGMD_MUTATION"]
    res["gene_name"] = gene_name
    res["ensembl_gene_id"] = ens_gene_id

    res = res[["ensembl_gene_id", "gene_name", "chromosome", "ref_allele", "start_pos", "alt_allele",
               "phenotype"]]
    res = res[res["phenotype"].notnull()]

    return res


def get_diseases_from_gene_name(gene_name, with_vars=False):
    """Retrieve diseases associated to a specific gene, using Ensembl.
    :param gene_name: [str] name of the query gene
    :param with_vars: [bool] True to also retrieve phenotypes associated to variants of the gene
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "location", "disease",
    "phenotypes"] or ["ensembl_gene_id", "gene_name", "location", "variation", "diseases",
    "phenotypes"] if with_vars=True
    """
    if gene_name.startswith("MT-"):
        gene_name = gene_name.upper().split("-")[1]

    server = "https://rest.ensembl.org"
    ext = "/phenotype/gene/homo_sapiens/{}?include_associated={}".format(gene_name, int(with_vars))
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    # if not r.ok:
    #     r.raise_for_status()
    #     print("Wrong request.")
    #     sys.exit()
    res = r.json()

    if with_vars:
        df = pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "location", "variation", "disease", "phenotypes"])
        try:
            ensembl_gene_id = Mitocarta.query.filter(Mitocarta.gene_symbol == gene_name).first().ensembl_id
        except AttributeError:  # gene not in Mitocarta
            return pd.DataFrame()
    else:
        df = pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "location", "disease", "phenotypes"])

    for el in res:
        if el == "error":
            return df

        if with_vars:
            if "attributes" in el and "associated_gene" in el["attributes"].keys() and "Variation" in el and "ontology_accessions" in el:
                row = pd.DataFrame({"ensembl_gene_id": ensembl_gene_id,
                                    "gene_name": [el["attributes"]["associated_gene"]],
                                    "location": el["location"],
                                    "variation": [el["Variation"]],
                                    "disease": [el["description"]],
                                    "phenotypes": [el["ontology_accessions"]]})
        else:
            if "ontology_accessions" in el:
                row = pd.DataFrame({"ensembl_gene_id": [el["Gene"]], "gene_name": gene_name,
                                    "location": el["location"], "disease": [el["description"]],
                                    "phenotypes": [el["ontology_accessions"]]})
        try:
            df = df.append(row, ignore_index=True)
        except UnboundLocalError:
            return df
    df.drop_duplicates(subset="disease", inplace=True)

    return df


def get_diseases_from_gene_id(ens_gene_id, with_vars=False):
    """
    Retrieve diseases associated with a specific Ensembl gene id, using Ensembl.
    :param ens_gene_id: [str] query Ensembl gene id
    :param with_vars: [bool] True to also retrieve phenotypes associated to variants of the gene
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "location", "disease",
    "phenotypes"]
    """
    try:
        gene_name = Mitocarta.query.filter(Mitocarta.ensembl_id == ens_gene_id).first().gene_symbol
    except AttributeError:  # gene not in Mitocarta
        return pd.DataFrame()

    # server = "https://rest.ensembl.org"
    # ext = "/phenotype/gene/homo_sapiens/{}?include_associated={}".format(gene_name, int(with_vars))
    # r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    # res = r.json()
    #
    # df = pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "location", "disease", "phenotypes"])
    # for el in res:
    #     row = pd.DataFrame({"ensembl_gene_id": [el["Gene"]], "gene_name": gene_name,
    #                         "location": el["location"], "disease": [el["description"]],
    #                         "phenotypes": [el["ontology_accessions"]]})
    #     df = df.append(row, ignore_index=True)
    # df.drop_duplicates(subset="disease", inplace=True)

    df = get_diseases_from_gene_name(gene_name, with_vars)

    return df


def final_from_gene_name(vars_df, disease_df, with_vars=False):
    """
    Create the final dataframe for gene data.
    :param vars_df: result of get_vars_from_gene_name()
    :param disease_df: result of get_diseases_from_gene_name()
    :return: pd.DataFrame with columns ["variant", "chromosome", "ensembl_gene_id", "gene_name",
    "variation", "disease_id", "disease", "phenotype_id", "phenotype_name"]
    """
    disease_df = disease_df[disease_df["location"].str.startswith(vars_df["chromosome"].value_counts().idxmax())]
    # TEST: following line might be useless: location refers to gene start and stop positions, not variants
    disease_df["position"] = disease_df["location"].str.split("-", expand=True)[0]
    vars_df["position"] = vars_df["chromosome"] + ":" + vars_df["start_pos"].astype(str)
    vars_df.drop(["ensembl_gene_id", "gene_name"], axis=1, inplace=True)
    df = (disease_df.set_index("position")
          .join(vars_df.set_index("position"))
          .reset_index())
    df["variant"] = df["ref_allele"] + df["start_pos"].astype(str) + df["alt_allele"]

    if with_vars:
        final_df = pd.DataFrame(columns=["variant", "chromosome", "ensembl_gene_id", "gene_name",
                                         "variation", "disease_id", "disease", "phenotype_id",
                                         "phenotype_name"])
    else:
        final_df = pd.DataFrame(columns=["variant", "chromosome", "ensembl_gene_id", "gene_name",
                                         "disease_id", "disease", "phenotype_id", "phenotype_name"])

    for row in df.itertuples():
        disease_id = disease_name_to_id(row.disease)
        for pheno in row.phenotypes:
            if with_vars:
                new_row = pd.DataFrame({"variant": row.variant, "chromosome": row.chromosome,
                                        "ensembl_gene_id": row.ensembl_gene_id,
                                        "gene_name": row.gene_name, "variation": row.variation,
                                        "disease_id": disease_id, "disease": row.disease,
                                        "phenotype_id": pheno,
                                        "phenotype_name": pheno_id_to_term(pheno)}, index=[row.variant])
            else:
                new_row = pd.DataFrame({"variant": row.variant, "chromosome": row.chromosome,
                                        "ensembl_gene_id": row.ensembl_gene_id,
                                        "gene_name": row.gene_name, "disease_id": disease_id,
                                        "disease": row.disease, "phenotype_id": pheno,
                                        "phenotype_name": pheno_id_to_term(pheno)},
                                       index=[row.variant])
            final_df = final_df.append(new_row, ignore_index=True)

    return final_df


def json_from_gene(gene_input):
    """
    Create the final json structure from gene data.
    :param gene_input: gene name to use for the queries
    :return:
    """
    vars_df = get_vars_from_gene_name(gene_input)
    disease_df = get_diseases_from_gene_name(gene_input, True)
    if disease_df.shape[0] == 0:
        disease_df = get_diseases_from_gene_name(gene_input)

    var_to_disease = VarDiseaseAss.query.filter(VarDiseaseAss.dbsnp_id.in_(vars_df.dbsnp_id)).all()
    for idx in vars_df.index:
        for item in var_to_disease:
            if fuzz.token_sort_ratio(vars_df.at[idx, "phenotype"].capitalize(),
                                     item.disease_name.capitalize()) >= 90:
                vars_df.at[idx, "phenotype"] = item.disease_name

    vars_df.drop(["ensembl_gene_id", "chromosome", "ref_allele", "start_pos",
                  "alt_allele"], axis=1, inplace=True)
    vars_df.rename({"phenotype": "disease_name"}, axis=1, inplace=True)
    vars_json = json.loads(vars_df.to_json(orient="records"))

    disease_df.drop(["ensembl_gene_id", "location"], axis=1, inplace=True)
    disease_df.rename({"disease": "disease_name", "phenotypes": "phenotype_ids"},
                      axis=1, inplace=True)

    gene_to_disease = GeneDiseaseAss.query.filter(GeneDiseaseAss.gene_symbol == gene_input).all()
    # TODO: it is possible to add more gene-related diseases from this table (but need to find phenotypes)
    for idx in disease_df.index:
        for item in gene_to_disease:
            if fuzz.token_sort_ratio(disease_df.at[idx, "disease_name"].capitalize(),
                                     item.disease_name.capitalize()) >= 90:
                disease_df.at[idx, "disease_name"] = item.disease_name

    disease_json = json.loads(disease_df.to_json(orient="records"))
    for el in disease_json:
        el["phenotype_names"] = []
        for pheno in el["phenotype_ids"]:
            el["phenotype_names"].append(pheno_id_to_term(pheno))
        el["disease_id"] = disease_name_to_id(el["disease_name"])  # TODO: correct this to grep from db table

    final_json = {}
    final_json["variants"] = vars_json
    final_json["diseases"] = disease_json

    return final_json


def network_from_gene_json(final_json):
    """
    Create the required nodes and edges dictionaries to build the network from gene data.
    :param final_json: output from json_from_gene()
    :return: dictionary with nodes and edges for network construction
    """
    var_json = final_json["variants"]
    dis_json = final_json["diseases"]

    v_variants = []
    v_variants.extend([el["variant"] for el in var_json])
    v_variants = set(v_variants)

    genes = []
    genes.extend([el["gene_name"] for el in var_json])
    genes = set(genes)

    v_diseases = []
    v_diseases.extend([el["disease_name"] for el in var_json])
    v_diseases = set(v_diseases)

    d_diseases = []
    d_diseases.extend([el["disease_name"] for el in dis_json])
    d_diseases = set(d_diseases)

    d_phenotypes = []
    for el in dis_json:
        d_phenotypes.extend(el["phenotype_names"])
    d_phenotypes = set(d_phenotypes)

    nodes = []
    edges = []
    ids = 0
    id_dict = {}
    for el in genes:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#739E82",
                                                        "border": "#5F826B"}})
        id_dict[el] = ids
    for el in v_variants:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#F9CF45",
                                                        "border": "#CCAA39"}})
        id_dict[el] = ids
    for el in v_diseases:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                        "border": "#B06A57"}})
        id_dict[el] = ids
    for el in d_diseases:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                        "border": "#B06A57"}})
        id_dict[el] = ids
    for el in d_phenotypes:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#93B5C6",
                                                        "border": "#7995A3"}})
        id_dict[el] = ids

    for el in var_json:
        # gene to variants
        edges.append({"from": id_dict[el["gene_name"]], "to": id_dict[el["variant"]]})
        # variant to diseases
        edges.append({"from": id_dict[el["variant"]], "to": id_dict[el["disease_name"]]})
    for el in dis_json:
        # gene to diseases
        edges.append({"from": id_dict[el["gene_name"]], "to": id_dict[el["disease_name"]]})
        # disease to phenotypes
        for pheno in el["phenotype_names"]:
            edges.append({"from": id_dict[el["disease_name"]], "to": id_dict[pheno]})

    return {"nodes": nodes, "edges": edges}


def network_from_gene(gene, vars_df, diseases_df):
    nodes = []
    edges = []
    ids = 1
    id_dict = {}
    # gene node
    nodes.append({"id": ids, "label": gene, "color": {"background": "#739E82",
                                                      "border": "#5F826B"}})
    id_dict[gene] = ids
    # disease nodes
    ass_diseases = GeneDiseaseAss.query.filter(GeneDiseaseAss.gene_symbol == gene,
                                               GeneDiseaseAss.score > 0.3).all()
    diseases = diseases_df.disease.unique().tolist()
    diseases.extend([el.disease_name for el in ass_diseases])
    for el in set(diseases):
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                        "border": "#B06A57"}})
        id_dict[el] = ids
    # variant nodes
    vars_df["variant"] = vars_df["ref_allele"] + vars_df["start_pos"].astype(str) + vars_df["alt_allele"]
    vars_df.drop_duplicates(subset=["variant", "phenotype"], inplace=True)
    for el in set(vars_df.variant):
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#F9CF45",
                                                        "border": "#CCAA39"}})
        id_dict[el] = ids
    # phenotype nodes
    phenos = []
    phenos_ids = []
    for el in diseases_df.phenotypes:
        phenos_ids.extend(el)
    for el in phenos_ids:
        phenos.append(pheno_id_to_term(el))
    phenos.extend(vars_df.phenotype.unique().tolist())
    for el in set(phenos):
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#93B5C6",
                                                        "border": "#7995A3"}})
        id_dict[el] = ids

    for el in vars_df.itertuples():
        # gene to variant edge
        edges.append({"from": id_dict[gene], "to": id_dict[el.variant]})
        # variant to phenotype edge
        edges.append({"from": id_dict[el.variant], "to": id_dict[el.phenotype]})
    for el in diseases_df.itertuples():
        # gene to disease edge
        edges.append({"from": id_dict[gene], "to": id_dict[el.disease]})
        # disease to phenotype edge
        for ph in el.phenotypes:
            edges.append({"from": id_dict[el.disease], "to": id_dict[pheno_id_to_term(ph)]})
        # edges.append({"from": id_dict[el.disease], "to": id_dict[el.]}) # TODO

    return {"nodes": nodes, "edges": edges}


def network_from_gene_name(final_df):
    """
    Create nodes and edges lists for the creation of the network, starting from a gene name.
    :param final_df: final dataframe returned by final_from_gene_name()
    :return: dictionary with a nodes list and an edges list
    """
    nodes = []
    edges = []
    ids = 0
    id_dict = {}
    genes = final_df["gene_name"].unique()
    for el in genes:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#739E82",
                                                        "border": "#5F826B"}})
        id_dict[el] = ids
    variants = final_df["variant"].unique()
    for el in variants:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#F9CF45",
                                                        "border": "#CCAA39"}})
        id_dict[el] = ids
    diseases = final_df["disease"].unique()
    for el in diseases:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#D7816A",
                                                        "border": "#B06A57"}})
        id_dict[el] = ids
    phenotypes = final_df["phenotype_name"].unique()
    for el in phenotypes:
        ids += 1
        nodes.append({"id": ids, "label": el, "color": {"background": "#93B5C6",
                                                        "border": "#7995A3"}})
        id_dict[el] = ids

    gene = genes[0]
    var_list = final_df[final_df["gene_name"] == gene].variant.unique()
    for dis in diseases:
        edges.append({"from": id_dict[gene], "to": id_dict[dis]})
        pheno_names = final_df[final_df["disease"] == dis].phenotype_name.unique()
        for pheno in pheno_names:
            edges.append({"from": id_dict[dis], "to": id_dict[pheno]})
    for var in var_list:
        edges.append({"from": id_dict[gene], "to": id_dict[var]})
        dis_names = final_df[final_df["variant"] == var].disease.unique()
        for dis in dis_names:
            edges.append({"from": id_dict[var], "to": id_dict[dis]})
            pheno_names = final_df[final_df["disease"] == dis].phenotype_name.unique()
            for pheno in pheno_names:
                edges.append({"from": id_dict[dis], "to": id_dict[pheno]})

    return {"nodes": nodes, "edges": edges}


# FROM PHENOTYPE #


def get_genes_from_phenotype(phenotype):
    """
    Retrieve genes related to a phenotype, using Ensembl.
    :param phenotype: [str] accession id of the phenotype to search for
    :return: pd.DataFrame with columns ["gene_name", "variation", "phenotypes", "description"]
    """
    server = "https://rest.ensembl.org"
    ext = "/phenotype/accession/homo_sapiens/{}?".format(phenotype)
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    # if not r.ok:
    #     r.raise_for_status()
    #     print("Wrong request.")
    #     sys.exit()
    res = r.json()

    # df = pd.DataFrame(columns=["gene_name", "variation", "phenotypes", "description"])
    df = pd.DataFrame(columns=["gene_name", "phenotypes", "description"])
    if len(res) != 0:
        for el in res:
            # if "attributes" in el and "associated_gene" in el["attributes"].keys() and "Variation" in el:
            if "attributes" in el and "associated_gene" in el["attributes"].keys():
                row = pd.DataFrame({"gene_name": [el["attributes"]["associated_gene"]],
                                    # "variation": [el["Variation"]],
                                    "phenotypes": [el["mapped_to_accession"]],
                                    "description": [el["description"]]})
            elif "attributes" in el and "Gene" in el:
                row = pd.DataFrame({"gene_name": [el["Gene"]],
                                    "phenotypes": [el["mapped_to_accession"]],
                                    "description": [el["description"]]})
            else:
                row = pd.DataFrame({"gene_name": [""],
                                    "phenotypes": [""],
                                    "description": [""]})
            df = df.append(row, ignore_index=True)
    else:
        dis_maps = DiseaseMappings.query.filter(DiseaseMappings.disease_id == phenotype).first()
        pheno_umls = dis_maps.umls_disease_id
        rel_genes = GeneDiseaseAss.query.filter(GeneDiseaseAss.umls_disease_id == pheno_umls).all()
        for el in rel_genes:
            row = pd.DataFrame({"gene_name": [el.gene_symbol], "phenotypes": [el.phenotype],
                                "description": [el.disease_name]})
            df = df.append(row, ignore_index=True)

    df.drop_duplicates(inplace=True)

    # TODO: find a better way to retrieve the list of Mitocarta genes from the table
    mito_genes = pd.read_sql("select * from Mitocarta", "sqlite:///hmtphenome.db")
    try:
        df["gene_name"] = df["gene_name"].str.split("-", expand=True)[1]
    except KeyError:
        pass
    df = df[df["gene_name"].isin(mito_genes["ensembl_id"])]
    # TODO: this check needs to be moved after the first part of the previous if, then falling back to the second part

    return df


def get_vars_from_phenotype(phenotype):
    """
    Retrieve variants related to a phenotype, exploiting the get_genes_from_phenotype() and
    get_vars_from_gene_name() functions.
    :param phenotype: [str] accession id of the phenotype to search for
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "chromosome", "ref_allele",
    "start_pos", "alt_allele", "phenotype", "phenotype_id"]
    """
    rel_genes = get_genes_from_phenotype(phenotype)
    if rel_genes.shape[0] == 0:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "chromosome", "ref_allele",
                                     "start_pos", "alt_allele", "phenotype", "phenotype_id"])
    try:
        pheno_names = rel_genes["description"].unique().tolist()
    except KeyError:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "chromosome", "ref_allele",
                                     "start_pos", "alt_allele", "phenotype", "phenotype_id"])

    rel_vars = pd.DataFrame()
    for el in set(rel_genes["gene_name"]):
        rel_vars = rel_vars.append(get_vars_from_gene_name(el))

    try:
        rel_vars = rel_vars[rel_vars["phenotype"].isin(pheno_names)]
        rel_vars["phenotype_id"] = phenotype
    except KeyError:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "chromosome", "ref_allele",
                                     "start_pos", "alt_allele", "phenotype", "phenotype_id"])

    return rel_vars


def get_diseases_from_phenotype(phenotype):
    """
    Retrieve diseases related to a phenotype, exploiting the get_genes_from_phenotype() and
    get_diseases_from_gene_name() functions.
    :param phenotype: [str] accession id of the phenotype to search for
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "location", "disease",
    "phenotype_ids"]
    """
    rel_genes = get_genes_from_phenotype(phenotype)
    # try:
    #     pheno_names = rel_genes["description"].unique().tolist()
    # except KeyError:
    #     return pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "chromosome", "ref_allele",
    #                                  "start_pos", "alt_allele", "phenotype"])

    rel_diseases = pd.DataFrame()
    for el in set(rel_genes["gene_name"]):
        rel_diseases = rel_diseases.append(get_diseases_from_gene_name(el))

    try:
        rel_diseases = rel_diseases[rel_diseases["phenotypes"].astype(str).str.contains(phenotype)]
        rel_diseases.rename({"phenotypes": "phenotype_ids"}, axis=1, inplace=True)
    except KeyError:
        return pd.DataFrame(columns=["ensembl_gene_id", "gene_name", "location", "disease",
                                     "phenotype_ids"])

    return rel_diseases


# TODO: I'm leaving this one as the last step, it needs more work
def final_from_phenotype(vars_df, disease_df):
    """
    Create the final dataframe for phenotype data.
    :param vars_df: result of get_vars_from_phenotype()
    :param disease_df: result of get_diseases_from_phenotype()
    :return: pd.DataFrame with columns ["variant", "chromosome", "ensembl_gene_id", "gene_name",
    "variation", "disease_id", "disease", "phenotype_id", "phenotype_name"]
    """
    disease_df = disease_df[disease_df["location"].str.startswith(vars_df["chromosome"].value_counts().idxmax())]
    disease_df["position"] = disease_df["location"].str.split("-", expand=True)[0]
    vars_df["position"] = vars_df["chromosome"] + ":" + vars_df["start_pos"].astype(str)
    vars_df.drop(["ensembl_gene_id", "gene_name"], axis=1, inplace=True)
    df = (disease_df.set_index("position")
          .join(vars_df.set_index("position"))
          .reset_index())
    df["variant"] = df["ref_allele"] + df["start_pos"].astype(str) + df["alt_allele"]

    final_df = pd.DataFrame(columns=["variant", "chromosome", "ensembl_gene_id", "gene_name",
                                     "variation", "disease_id", "disease", "phenotype_id",
                                     "phenotype_name"])

    for row in df.itertuples():
        disease_id = disease_name_to_id(row.disease)
        for pheno in row.phenotypes:
            new_row = pd.DataFrame({"variant": row.variant, "chromosome": row.chromosome,
                                    "ensembl_gene_id": row.ensembl_gene_id,
                                    "gene_name": row.gene_name, "variation": row.variation,
                                    "disease_id": disease_id, "disease": row.disease,
                                    "phenotype_id": pheno,
                                    "phenotype_name": pheno_id_to_term(pheno)}, index=[row.variant])
            final_df = final_df.append(new_row, ignore_index=True)

    return final_df


# FROM DISEASE #


def get_genes_from_disease_name(disease_name):
    """
    Retrieve genes related to a disease, using Biomart.
    :param disease_name: [str] name of the query disease
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "gene_descr", "chromosome",
    "start_pos", "stop_pos", "disease"]
    """
    server = Server(host="http://www.ensembl.org")
    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets["hsapiens_gene_ensembl"]
    res = dataset.query(attributes=["ensembl_gene_id", "external_gene_name", "description",
                                    "chromosome_name", "start_position", "end_position",
                                    "phenotype_description"],
                        filters={"phenotype_description": disease_name})

    res.rename({"Gene stable ID": "ensembl_gene_id", "Gene name": "gene_name",
                "Gene description": "gene_descr", "Chromosome/scaffold name": "chromosome",
                "Gene start (bp)": "start_pos", "Gene end (bp)": "stop_pos",
                "Phenotype description": "disease"}, axis=1, inplace=True)
    res.drop_duplicates(inplace=True)

    return res


def get_vars_from_disease_name(disease_name):
    """
    Retrieve variants related to a specific disease, exploiting the get_genes_from_disease_name()
    and get_vars_from_gene_name() functions.
    :param disease_name: [str] name of the query disease
    :return: pd.DataFrame with columns ["gene", "ensembl_gene_id", "chromosome", "ref_allele",
    "start_pos", "alt_allele", "disease"]
    """
    rel_genes = get_genes_from_disease_name(disease_name)
    try:
        gene_ids = rel_genes["ensembl_gene_id"].unique()
    except KeyError:
        return pd.DataFrame()

    rel_vars = pd.DataFrame()
    for el in set(gene_ids):
        rel_vars = rel_vars.append(get_vars_from_gene_id(el))

    try:
        rel_vars = rel_vars[rel_vars["ensembl_gene_id"].isin(gene_ids)]
    except KeyError:
        return pd.DataFrame(columns=["gene", "ensembl_gene_id", "chromosome", "ref_allele",
                                     "start_pos", "alt_allele", "disease"])

    # TODO: find a better way to retrieve variants related to a specific disease
    rel_vars = rel_vars[rel_vars["phenotype"] == disease_name]
    rel_vars.rename({"phenotype": "disease"}, axis=1, inplace=True)

    return rel_vars








