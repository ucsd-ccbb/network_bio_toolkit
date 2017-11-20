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
                                './slowkow_databases/Marbach2016_TF.txt'], sep = '\n'):

    """
        This function assumes you are using files that are formatted as a single list of "sep" separated
        transcription factors. While this method is intended to load the slowkow database, any files of the
        correct format can be processed using this method.

        Args:
            filename_list: list of strings, a list of six files from which to extract the slowkow database, or whichever files you wish
                to extract transcription factors from
            sep: string, the type of deliminator separating your TF entries

        Returns: A list of all unique transcription factors from the lists provided

    """

    # read files formatted as \n separated items
    return_list = []
    for file_name in filename_list:
        with open(file_name) as f:
            lines = f.read().split(sep)
            return_list.extend(lines)

    # convert everything to ALL CAPS
    [x.upper() for x in return_list]

    # remove duplicates
    return set(return_list)




def load_jaspar(filename = 'jaspar_genereg_matrix.txt'):

    """
        This function loads The 2016 version of Jaspar's TF database, which can be dowloaded as a .txt file from
        http://jaspar2016.genereg.net/. At the bottom of the page, navigate to the "Download" button. Open the
        "database" file and download "MATRIX.txt".

        Args:
            filename = string, the path to the jaspar "MATRIX.txt" file

        Returns: A list of all unique transcription factors from the file provided

    """

    # parse jaspar file
    jasp_df = pd.read_csv(filename, sep="\t", header=None, names=['col1', 'col2', 'col3', 'col4', 'tf_genes'])

    # return transcription factors with ALL CAPS names
    return list(set(jasp_df['tf_genes'].str.upper()))



def create_TF_list(slowkow_bool=True,
                    slowkow_files=['./slowkow_databases/TRED_TF.txt',
                                      './slowkow_databases/ITFP_TF.txt',
                                      './slowkow_databases/ENCODE_TF.txt',
                                      './slowkow_databases/Neph2012_TF.txt',
                                      './slowkow_databases/TRRUST_TF.txt',
                                      './slowkow_databases/Marbach2016_TF.txt'],
                    slowkow_sep = '\n',
                    jaspar_bool=True,
                    jaspar_file="jaspar_genereg_matrix.txt"):
    """
        This function loads The 2016 version of Jaspar's TF database, which can be dowloaded as a .txt file from
        http://jaspar2016.genereg.net/. At the bottom of the page, navigate to the "Download" button. Open the
        "database" file and download "MATRIX.txt".

        Args:
            slowkow_bool = boolean, whether or not to include the slowkow databases
            slowkow_files = list or strings, a list of input file path names to the load_slowkow method
            slowkow_sep = string, the type of deliminator used in the slowkow_files to separate the transcription factors
            jaspar_bool = boolean, whether or not to include the jaspar database
            jaspar_file = string, the file path where to find the jaspar database

        Returns: A list of all unique transcription factors from the files provided

    """

    TF_list = []

    if slowkow_bool == True:
        slowkow_TFs = load_slowkow(slowkow_files, slowkow_sep)
        TF_list.extend(slowkow_TFs)

    if jaspar_bool == True:
        jaspar_TFs = load_jaspar(jaspar_file)
        TF_list.extend(jaspar_TFs)

    return list(set(TF_list))


def load_and_process_small_STRING(filename="STRING_network.xlsx"):

    """
        This function loads a version of the STRING database to serve as our background network. It converts the STRING
        excel file into a pandas dataframe and extracts the edge information from that dataframe. We create two dictionaries
        to be used when creating a networkx graph: one maps edges to their weights, and the other maps edges to their
        weight's sign (+ or -)

        Args:
            filename = string, the path to the STRING excel file

        Returns: STRING_DF = Dataframe, a pandas dataframe of the entire STRING file
                 db_edges = list, a list of tuples of the form (source, target, edge weight)
                 db_sign_att = list, a list of tuples of the form (source, target, sign of edge weight)

    """

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

    """
        This function takes as input two lists (as created by our load_and_process_STRING function), and removes
        all of the edges whose source is not a TF designated in our TF list (created by a function such as create_TF_list).

        Args:
            db_edges = list, a list of tuples of the form (source, target, edge weight)
            db_sign_att = list, a list of tuples of the form (source, target, sign of edge weight)
            TF_list = list, a list of transcription factors produced from a function such as create_TF_list

        Returns: edge_list_filtered = list of same form as db_edges, where only TF-as-source edges remain
                 sign_att_list_filtered = list of same form as db_sign_att, where only TF-as-source edges remain

    """

    # extracting TR edge information from background database
    edge_list_filtered = []
    sign_att_list_filtered = []
    for i in range(len(db_edges)):
        if db_edges[i][0] in list(TF_list):
            edge_list_filtered.append(db_edges[i])
            sign_att_list_filtered.append(db_sign_att[i])

    return edge_list_filtered, sign_att_list_filtered


def make_digraph(db_edges, db_sign_att, TF_list):

    """
        This function filters down the to given lists of edges-and-attribute by calling filter_background. This removes
        all edge from our lists whose source node is not a transcription factor (as specified by TF_list). Our function
        then creates a networkX directed graph fom these two edges-and-attribute lists

        Args:
            db_edges = list, a list of tuples of the form (source, target, edge weight)
            db_sign_att = list, a list of tuples of the form (source, target, sign of edge weight)
            TF_list = list, a list of transcription factors produced from a function such as create_TF_list

        Returns: DG = DiGraph, a networkX directed graph with only TF-as-source nodes, containing attributes 'weight'
                      and 'sign'
    """

    # use only edges from background network associated with our TF list
    edge_list_filtered, sign_att_list_filtered = filter_background(db_edges, db_sign_att, TF_list)

    # create networkx digraph from weighted edge list, add sign edge attributes
    DG = nx.DiGraph()
    DG.add_weighted_edges_from(edge_list_filtered)
    for i in range(len(sign_att_list_filtered)):
        DG[sign_att_list_filtered[i][0]][sign_att_list_filtered[i][1]]['sign'] = sign_att_list_filtered[i][2]

    return DG


def load_DEG_with_up_downs(filename="differencially_expressed_genes.txt", filter_value=0.3):

    """
        This function loads up a list of differentially expressed genes, removes genes with lfdr value greater
        than the specified filter value, then extracts the up/down regulation information associated with those genes.

        Args:
            filename = string, path to our list of differencially expressed genes, containing up/down regulated information
                       under column log2.89.12
            filter_value = float, cutoff lfdr value to filter our DEG's by

        Returns: DEG_list = list, a list of differencially expressed genes that have a lfdr value less than the filter value
                 DEG_to_updown = dict, a dictionary that maps a DEG name to wether it is up regulated or down regulated. If
                                 no up/down regulation info, will set value to zero
    """

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
                DEG_to_updown[str(DEG_db.symbol[i]).upper()] = DEG_db['log2.89.12'][i]

    return DEG_list, DEG_to_updown


def add_updown_from_DEG(DG, DEG_filename="differencially_expressed_genes.txt", DEG_filter_value=0.3):

    """
        This function loads up a list of differentially expressed genes, removes genes with lfdr value greater
        than the specified filter value, then extracts the up/down regulation information associated with those genes.
        Our function then finds all DEG's located withing graph DG, and adds these up/down regulation values as a node
        attribute.

        Args:
            DG = DiGraph, a networkX directed graph
            DEG_filename = string, path to our list of differencially expressed genes, containing up/down regulated information
                       under column log2.89.12
            DEG_filter_value = float, cutoff lfdr value to filter our DEG's by

        Returns: DEG_list = list, a list of differencially expressed genes that have a lfdr value less than the filter value
    """

    DEG_list, DEG_to_updown = load_DEG_with_up_downs(DEG_filename, DEG_filter_value)

    # get all the differencially expressed genes in DG
    DEG_in_DG = set(DG.nodes()) & set(DEG_list)

    # add node attribute to each node in DG if it exists, otherwise set to zero
    zero_dict = dict(zip(DG.nodes(), [0] * len(DG.nodes())))
    for gene in DEG_in_DG:
        zero_dict[gene] = DEG_to_updown[gene]
    nx.set_node_attributes(DG, 'updown', zero_dict)

    return DEG_list


