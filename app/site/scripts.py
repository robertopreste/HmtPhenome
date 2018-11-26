#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
from app.site.models import Mitocarta, Phenotypes, Diseases
import pandas as pd
import numpy as np
import requests
# import sys
from pybiomart import Server


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
        pheno_name = Phenotypes.query.filter(Phenotypes.hpo_id == pheno_id).first().hpo_term_name
    else:
        efo_id = pheno_id.strip("EFO:")
        base_url = "https://www.ebi.ac.uk/ols/api/ontologies/efo/terms?"
        iri_url = "iri=http://www.ebi.ac.uk/efo/EFO_{}".format(efo_id)
        r = requests.get(base_url + iri_url, headers={"Content-Type": "application/json"})
        res = r.json()
        pheno_name = res["_embedded"]["terms"][0]["label"]

    return pheno_name


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
    :return:
    """
    gene_name = get_gene_from_variant(chrom, var_start, var_end)["gene_name"][0]
    diseases = get_diseases_from_gene_name(gene_name, True)

    if var_end is not None:
        diseases = diseases[(diseases["location"].str.startswith(chrom)) &
                            (diseases["location"].str.contains(var_start)) &
                            (diseases["location"].str.contains(var_end))]
    else:
        diseases = diseases[(diseases["location"].str.startswith(chrom)) &
                            (diseases["location"].str.contains(var_start))]

    return diseases


def final_from_variant(gene_df, pheno_df, disease_df):
    """
    Create the final dataframe for variant data.
    :param gene_df: result of get_gene_from_variant()
    :param pheno_df: result of get_pheno_from_variant()
    :param disease_df: result of get_diseases_from_variant()
    :return: pd.DataFrame with columns ["variant", "ensembl_gene_id", "gene_name", "variation",
    "disease", "phenotype_id", "phenotype_name"]
    """
    df = (disease_df.set_index("disease")
          .join(pheno_df.set_index("phenotype"))
          .reset_index())
    df["variant"] = df["ref_allele"] + df["start_pos"].astype(str) + df["alt_allele"]
    # full_df = full_df[["ensembl_gene_id", "gene_name", "variation", "variant", "disease",
    #                    "phenotypes"]]
    final_df = pd.DataFrame(columns=["variant", "ensembl_gene_id", "gene_name", "variation",
                                     "disease", "phenotype_id", "phenotype_name"])

    for row in df.itertuples():
        for pheno in row.phenotypes:
            new_row = pd.DataFrame({"variant": row.variant, "ensembl_gene_id": row.ensembl_gene_id,
                                    "gene_name": row.gene_name, "variation": row.variation,
                                    "disease": row.disease, "phenotype_id": pheno,
                                    "phenotype_name": pheno_id_to_term(pheno)}, index=[row.variant])
            final_df = final_df.append(new_row, ignore_index=True)

    return final_df


def get_vars_from_gene_name(gene_name):
    """
    Retrieve all variants associated with a specific gene, using Biomart.
    :param gene_name: [str] name of the query gene
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "chromosome", "ref_allele",
    "start_pos", "alt_allele", "phenotype"]
    """
    gene_name = gene_name.upper().lstrip("MT-")

    try:
        ens_gene_id = Mitocarta.query.filter(Mitocarta.gene_symbol == gene_name).first().ensembl_id
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

    return res


def get_diseases_from_gene_name(gene_name, with_vars=False):
    """Retrieve diseases associated to a specific gene, using Ensembl.
    :param gene_name: [str] name of the query gene
    :param with_vars: [bool] True to also retrieve phenotypes associated to variants of the gene
    :return: pd.DataFrame with columns ["ensembl_gene_id", "gene_name", "location", "disease",
    "phenotypes"] or ["ensembl_gene_id", "gene_name", "location", "variation", "diseases",
    "phenotypes"] if with_vars=True
    """
    gene_name = gene_name.upper().lstrip("MT-")

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

    df = pd.DataFrame(columns=["gene_name", "variation", "phenotypes", "description"])
    for el in res:
        if "attributes" in el and "associated_gene" in el["attributes"].keys() and "Variation" in el:
            row = pd.DataFrame({"gene_name": [el["attributes"]["associated_gene"]],
                                "variation": [el["Variation"]],
                                "phenotypes": [el["mapped_to_accession"]],
                                "description": [el["description"]]})
            df = df.append(row, ignore_index=True)
    df.drop_duplicates(inplace=True)

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


def disease_id_to_name(disease_id):
    """
    Convert a given disease ID into its common name.
    :param disease_id: [str] query disease ID
    :return: related disease name
    """
    q = Diseases.query.filter(Diseases.disease_id == disease_id).first()

    return q.disease_name


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




