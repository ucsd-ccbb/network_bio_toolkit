"""
-------------------------------------------
Author: Mikayla Webster (m1webste@ucsd.edu)
Date: 10/13/17
-------------------------------------------
"""

import pandas as pd
import networkx as nx


def load_slowkow(filename_list=['./slowkow_databases/TRED_TF.txt',
                                './slowkow_databases/ITFP_TF.txt',
                                './slowkow_databases/ENCODE_TF.txt',
                                './slowkow_databases/Neph2012_TF.txt',
                                './slowkow_databases/TRRUST_TF.txt',
                                './slowkow_databases/Marbach2016_TF.txt']):


    # read files formatted as \n separated items
    return_list = []
    for file_name in filename_list:
        with open(file_name) as f:
            lines = f.read().splitlines()
            return_list.extend(lines)

    # convert everything to ALL CAPS
    [x.upper() for x in return_list]

    # remove duplicates
    return set(return_list)




def load_jaspar(filename):
    # parse jaspar file
    jasp_df = pd.read_csv(filename, sep="\t", header=None, names=['col1', 'col2', 'col3', 'col4', 'tf_genes'])

    # return transcription factors with ALL CAPS names
    return list(jasp_df['tf_genes'].str.upper())



def create_TF_list(slowkow_bool=True,
                       slowkow_files=['./slowkow_databases/TRED_TF.txt',
                                      './slowkow_databases/ITFP_TF.txt',
                                      './slowkow_databases/ENCODE_TF.txt',
                                      './slowkow_databases/Neph2012_TF.txt',
                                      './slowkow_databases/TRRUST_TF.txt',
                                      './slowkow_databases/Marbach2016_TF.txt'],
                       jaspar_bool=True,
                       jaspar_file="jaspar_genereg_matrix.txt"):


    TF_list = []

    if slowkow_bool == True:
        slowkow_TFs = load_slowkow(slowkow_files)
        TF_list.extend(slowkow_TFs)

    if jaspar_bool == True:
        jaspar_TFs = load_jaspar(jaspar_file)
        TF_list.extend(jaspar_TFs)

    return list(set(TF_list))


def load_and_process_STRING(filename="STRING_network.xlsx"):
    # Load STRING database as background network
    STRING_DF = pd.read_excel(filename)

    # convert sources (sym1) and targets (sym2) to all caps
    sym1_list = STRING_DF.Symbol1.str.upper()
    sym2_list = STRING_DF.Symbol2.str.upper()

    # make an edge list with associated edge weight (db_edges)
    weight_list = STRING_DF.Weight
    db_edges = zip(sym1_list, sym2_list, weight_list)

    # make an edge list with associated activating (+)/inhibiting (-) sign (db_sign_att)
    sign_list = STRING_DF.Edge_Sign
    sign_num_list = []
    for sign in sign_list:
        if str(sign) == '+':
            sign_num_list.append(1)
        elif str(sign) == '-':
            sign_num_list.append(-1)
        else:
            sign_num_list.append(0)
    db_sign_att = zip(sym1_list, sym2_list, sign_num_list)

    return STRING_DF, db_edges, db_sign_att


def filter_background(db_edges, db_sign_att, TF_list):
    # extracting TR edge information from background database
    edge_list_filtered = []
    sign_att_list_filtered = []
    for i in range(len(db_edges)):
        if db_edges[i][0] in list(TF_list):
            edge_list_filtered.append(db_edges[i])
            sign_att_list_filtered.append(db_sign_att[i])

    return edge_list_filtered, sign_att_list_filtered


def make_digraph(db_edges, db_sign_att, TF_list):
    # use only edges from background network associated with our TF list
    edge_list_filtered, sign_att_list_filtered = filter_background(db_edges, db_sign_att, TF_list)

    # create networkx digraph from weighted edge list, add sign edge attributes
    DG = nx.DiGraph()
    DG.add_weighted_edges_from(edge_list_filtered)
    for i in range(len(sign_att_list_filtered)):
        DG[sign_att_list_filtered[i][0]][sign_att_list_filtered[i][1]]['sign'] = sign_att_list_filtered[i][2]

    return DG


def load_DEG_with_up_downs(filename="differencially_expressed_genes.txt", filter_value=0.3):
    # load differencially expressed genes (experimental results)
    DEG_db = pd.read_csv(filename, sep="\t")

    # filtering for lfdr < 0.3
    DEG_list = []
    DEG_to_updown = {}
    for i in range(len(DEG_db)):

        # removing Nan values
        if str(DEG_db.symbol[i]).upper() != 'NAN':

            # filtering DEG list by lfdr < filter_value
            if (DEG_db['lfdr.89.12'][i] < filter_value):
                DEG_list.append(str(DEG_db.symbol[i]).upper())

                # creating dictionary between DEG symbols and their up/down value
                if DEG_db['log2.89.12'][i] != 0:
                    DEG_to_updown[str(DEG_db.symbol[i]).upper()] = DEG_db['log2.89.12'][i]
                else:
                    DEG_to_updown[str(DEG_db.symbol[i]).upper()] = 0

    return DEG_list, DEG_to_updown


def add_updown_from_DEG(DG, DEG_filename="differencially_expressed_genes.txt", DEG_filter_value=0.3):
    DEG_list, DEG_to_updown = load_DEG_with_up_downs(DEG_filename, DEG_filter_value)

    # get all the differencially expressed genes in DG
    DEG_in_DG = set(DG.nodes()) & set(DEG_list)

    # add node attribute to each node in DG if it exists, otherwise set to zero
    zero_dict = dict(zip(DG.nodes(), [0] * len(DG.nodes())))
    for gene in DEG_in_DG:
        zero_dict[gene] = DEG_to_updown[gene]
    nx.set_node_attributes(DG, 'updown', zero_dict)

    return DEG_list