#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste
import os
import pandas as pd
import wget
import xml.etree.cElementTree as et


def get_node_val(node):
    """
    Return the value of the given node, if present, otherwise return None.
    :param node: node provided by cElementTree.getroot()
    :return:
    """
    return node.text if node is not None else None


def download_source(url, out):
    """
    Retrieve the given file from the web and save it with the given name in data/raw/.
    :param url: URL of the file to download
    :param out: name of the output file
    :return:
    """
    print("Downloading file from {} and saving it to {}...".format(url, out))
    wget.download(url, "data/raw/{}".format(out))
    print("Complete.")


def process_mitocarta(in_file, out_file):

    print("Processing Mitocarta data...")
    mitocarta = pd.read_excel("data/raw/{}".format(in_file), sheet_name=1)
    mitocarta = mitocarta[["EnsemblGeneID", "Symbol", "Description", "hg19_Chromosome",
                           "hg19_Start", "hg19_Stop"]]
    mitocarta.rename({"EnsemblGeneID": "ensembl_id", "Symbol": "gene_symbol",
                      "Description": "description", "hg19_Chromosome": "hg_chr",
                      "hg19_Start": "hg_start", "hg19_Stop": "hg_stop"}, axis=1, inplace=True)

    # Add mitochondrial tRNA genes to the list of Mitocarta genes
    mt_trnas = """ENSG00000210127,TA,tRNA Alanine,chrM,5587,5655
    ENSG00000210140,TC,tRNA Cysteine,chrM,5761,5826
    ENSG00000210154,TD,tRNA Aspartic Acid,chrM,7518,7585
    ENSG00000210194,TE,tRNA Glutamic Acid,chrM,14674,14742
    ENSG00000210049,TF,tRNA Phenylalanine,chrM,577,647
    ENSG00000210164,TG,tRNA Glycine,chrM,9991,10058
    ENSG00000210176,TH,tRNA Histidine,chrM,12138,12206
    ENSG00000210100,TI,tRNA Isoleucine,chrM,4263,4331
    ENSG00000210156,TK,tRNA Lysine,chrM,8295,8364
    ENSG00000209082,TL1,tRNA Leucine 1,chrM,3230,3304
    ENSG00000210191,TL2,tRNA Leucine 2,chrM,12266,12336
    ENSG00000210112,TM,tRNA Methionine,chrM,4402,4469
    ENSG00000210135,TN,tRNA Asparagine,chrM,5657,5729
    ENSG00000210196,TP,tRNA Proline,chrM,15956,16023
    ENSG00000210107,TQ,tRNA Glutamine,chrM,4329,4400
    ENSG00000210174,TR,tRNA Arginine,chrM,10405,10469
    ENSG00000210151,TS1,tRNA Serine 1,chrM,7446,7514
    ENSG00000210184,TS2,tRNA Serine 2,chrM,12207,12265
    ENSG00000210195,TT,tRNA Threonine,chrM,15888,15953
    ENSG00000210077,TV,tRNA Valine,chrM,1602,1670
    ENSG00000210117,TW,tRNA Tryptophan,chrM,5512,5579
    ENSG00000210144,TY,tRNA Tyrosine,chrM,5826,5891"""

    mt_trna_list = mt_trnas.split("\n")

    for el in mt_trna_list:
        row = el.split(",")
        new_row = pd.DataFrame({"ensembl_id": [row[0]], "gene_symbol": [row[1]],
                                "description": [row[2]], "hg_chr": [row[3]], "hg_start": [row[4]],
                                "hg_stop": [row[5]]})
        mitocarta = mitocarta.append(new_row, ignore_index=True)

    print("Saving processed Mitocarta data...")
    mitocarta.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def process_hpo_disgenphen(in_file, out_file):
    print("Processing HPO disease_gene_pheno data...")
    hpo = pd.read_csv("data/raw/{}".format(in_file), sep="\t", skiprows=1,
                      names=["disease_id", "gene_symbol", "entrez_gene_id", "hpo_id",
                             "hpo_term_name"])
    print("Saving processed HPO disease_gene_pheno data...")
    hpo.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def process_hpo_genphen(in_file, out_file):
    print("Processing HPO gene_pheno data...")
    hpo = pd.read_csv("data/raw/{}".format(in_file), sep="\t", skiprows=1,
                      names=["entrez_gene_id", "gene_symbol", "hpo_term_name", "hpo_id"])
    print("Saving processed HPO gene_pheno data...")
    hpo.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def process_hpo_phengen(in_file, out_file):
    print("Processing HPO pheno_gene data...")
    hpo = pd.read_csv("data/raw/{}".format(in_file), sep="\t", skiprows=1,
                      names=["hpo_id", "hpo_term_name", "entrez_gene_id", "gene_symbol"])
    print("Saving processed HPO pheno_gene data...")
    hpo.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def process_omim(in_file, out_file):
    print("Processing OMIM data...")
    omim = pd.read_csv("data/raw/{}".format(in_file), sep="\t", skiprows=3, comment="#",
                       names=["prefix", "mim_number", "mim_name", "alter_title", "includ_title"])
    omim = omim[["mim_number", "mim_name", "prefix"]]
    print("Saving processed OMIM data...")
    omim.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def process_orphanet(in_file, out_file):
    print("Processing Orphanet data...")
    parsed_xml = et.parse("data/raw/{}".format(in_file))
    orpha = pd.DataFrame(columns=["orpha_name", "orpha_num"])
    for node in parsed_xml.getroot():
        for el in node:
            dis_name = el.find("Name")
            orpha_num = el.find("OrphaNumber")
            orpha = orpha.append(pd.DataFrame({"orpha_name": [get_node_val(dis_name)],
                                               "orpha_num": [get_node_val(orpha_num)]}),
                                 ignore_index=True)
    print("Saving processed Orphanet data...")
    orpha.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def process_disgenet_gene(in_file, out_file, mitocarta_file):
    print("Processing Disgenet gene_disease data...")
    disge = pd.read_csv("data/raw/{}".format(in_file), sep="\t", encoding="latin1",
                        names=["entrez_gene_id", "gene_symbol", "umls_disease_id", "disease_name",
                               "score", "nofpmids", "nofsnps", "source"])
    disge = disge[["entrez_gene_id", "gene_symbol", "umls_disease_id", "disease_name", "score"]]
    mitocarta = pd.read_csv("data/tables/{}".format(mitocarta_file))
    disge = disge[disge.gene_symbol.isin(mitocarta.gene_symbol)]
    print("Saving processed Disgenet gene_disease data...")
    disge.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def process_disgenet_vars(in_file, out_file):
    print("Processing Disgenet vars_disease data...")
    disge = pd.read_csv("data/raw/{}".format(in_file), sep="\t", encoding="latin1",
                        names=["dbsnp_id", "umls_disease_id", "disease_name", "score", "nofpmids",
                               "source"])
    disge = disge[["dbsnp_id", "umls_disease_id", "disease_name", "score"]]
    print("Saving processed Disgenet vars_disease data...")
    disge.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def process_disgenet_maps(in_file, out_file):
    print("Processing Disgenet disease_mappings data...")
    disge = pd.read_csv("data/raw/{}".format(in_file), sep="\t",
                        names=["umls_disease_id", "disease_name", "vocabulary", "disease_id",
                               "alt_disease_name"])
    print("Saving processed Disgenet disease_mappings data...")
    disge.to_csv("data/tables/{}".format(out_file), index=False)
    print("Complete.")


def create_diseases(hpo_file, omim_file, orpha_file, out_file):
    print("Creating Diseases table...")
    hpo = pd.read_csv("data/tables/{}".format(hpo_file))
    omim = pd.read_csv("data/tables/{}".format(omim_file))
    orpha = pd.read_csv("data/tables/{}".format(orpha_file))

    dis_ids = hpo["disease_id"].unique()
    diseases = pd.DataFrame(columns=["disease_id", "disease_name"])
    for el in dis_ids:
        resource, resource_id = el.split(":")
        if resource == "OMIM":
            disease_name = omim[omim.mim_number == int(resource_id)]["mim_name"].values[0]
        elif resource == "ORPHA":
            disease_name = orpha[orpha.orpha_num == int(resource_id)]["orpha_name"].values[0]
        else:
            continue
        diseases = diseases.append(pd.DataFrame({"disease_id": [el], "disease_name": [disease_name]}))
    print("Saving Diseases table...")
    diseases.to_csv("data/tables/{}".format(out_file))
    print("Complete.")


def create_phenotypes(hpo_file, out_file):
    print("Creating Phenotypes table...")
    hpo = pd.read_csv("data/tables/{}".format(hpo_file))
    phenos = hpo[["hpo_id", "hpo_term_name"]]
    phenos.drop_duplicates(inplace=True)
    print("Saving Phenotypes table...")
    phenos.to_csv("data/tables/{}".format(out_file))
    print("Complete.")


if __name__ == '__main__':
    sources = {"mitocarta": ("http://www.broadinstitute.org/ftp/distribution/metabolic/papers/Pagliarini/MitoCarta2.0/Human.MitoCarta2.0.xls",
                             "mitocarta.xls",
                             "Mitocarta.csv"),
               "hpo_1": ("http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt",
                         "hpo_disease_gene_phenotype.txt",
                         "HpoDisGenePhen.csv"),
               "hpo_2": ("http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt",
                         "hpo_gene_phenotype.txt",
                         "HpoGenePhen.csv"),
               "hpo_3": ("http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastStableBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt",
                         "hpo_phenotype_gene.txt",
                         "HpoPhenGene.csv"),
               "omim": ("https://data.omim.org/downloads/IrNkLfaiTB6CrGnAGCXDzw/mimTitles.txt",
                        "omim_data.txt",  # TODO: change this OMIM API key because it won't work
                        "Omim.csv"),
               "orpha": ("http://www.orphadata.org/data/xml/en_product1.xml",
                         "orphanet_data.xml",
                         "Orphanet.csv"),
               "disge_gene": ("http://www.disgenet.org/ds/DisGeNET/results/all_gene_disease_associations.tsv.gz",
                              "disgenet_all_gene_disease.tsv.gz",
                              "GeneDiseaseAss.csv"),
               "disge_vars": ("http://www.disgenet.org/ds/DisGeNET/results/all_variant_disease_associations.tsv.gz",
                              "disgenet_all_vars_disease.tsv.gz",
                              "VarDiseaseAss.csv"),
               "disge_maps": ("http://www.disgenet.org/ds/DisGeNET/results/disease_mappings.tsv.gz",
                              "disgenet_mappings.tsv.gz",
                              "DiseaseMappings.csv")}

    if not os.path.isdir(os.path.join(os.getcwd(), "data/raw/")):
        os.makedirs(os.path.join(os.getcwd(), "data/raw/"))
    if not os.path.isdir(os.path.join(os.getcwd(), "data/tables/")):
        os.makedirs(os.path.join(os.getcwd(), "data/tables/"))

    for el in sources:
        download_source(sources[el][0], sources[el][1])

    process_mitocarta(sources["mitocarta"][1], sources["mitocarta"][2])
    process_hpo_disgenphen(sources["hpo_1"][1], sources["hpo_1"][2])
    process_hpo_genphen(sources["hpo_2"][1], sources["hpo_2"][2])
    process_hpo_phengen(sources["hpo_3"][1], sources["hpo_3"][2])
    process_omim(sources["omim"][1], sources["omim"][2])
    process_orphanet(sources["orpha"][1], sources["orpha"][2])
    process_disgenet_gene(sources["disge_gene"][1], sources["disge_gene"][2], sources["mitocarta"])
    process_disgenet_vars(sources["disge_vars"][1], sources["disge_vars"][2])
    process_disgenet_maps(sources["disge_maps"][1], sources["disge_maps"][2])
    create_diseases(sources["hpo_1"][2], sources["omim"][2], sources["orpha"][2], "Diseases.csv")
    create_phenotypes(sources["hpo_1"][2], "Phenotypes.csv")



