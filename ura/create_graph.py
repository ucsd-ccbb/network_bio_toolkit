"""
-------------------------------------------
Author: Mikayla Webster (m1webste@ucsd.edu)
Date: 10/13/17
-------------------------------------------
"""

import pandas as pd
import networkx as nx
import mygene


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

# Loss of 110/about 40000 edges do to removal of multiedges
def load_small_STRING_to_digraph(filename="../STRING_network.xlsx"):

    # Load STRING database as background network
    STRING_DF = pd.read_excel(filename)
    STRING_DF.drop(['ID1','ID2','Source_Type','Sign_Score','NumberOfScreens','Interaction_Database'], axis=1, inplace=True)
    STRING_DF.rename(index=str, columns={'Symbol1':'source', 'Symbol2':'target','Edge_Sign':'sign','Weight':'weight'}, inplace = True)

    # make all gene symbol names upper case
    STRING_DF = pd.concat([STRING_DF[col].astype(str).str.upper() for col in STRING_DF.columns], axis=1)

    # mark sign attribute as +1 for activating and -1 for inhibiting
    STRING_DF[STRING_DF == '+'] = 1
    STRING_DF[STRING_DF == '-'] = -1

    # make digraph
    STRING_DF['weight'] = map(lambda x: float(x), STRING_DF['weight'])
    G_str = nx.from_pandas_dataframe(STRING_DF, 'source', 'target', ['sign', 'weight'], create_using=nx.DiGraph())

    return G_str


def load_STRING_to_digraph(filename = "9606.protein.actions.v10.5.txt",confidence_filter=400):

    # read STRING file to dataframe
    df_full = pd.read_csv(filename, sep="\t")

    df = df_full.loc[df_full['score'] > confidence_filter]  # filter by confidence
    df_act = df.loc[df['action'] == 'activation']  # filter by activation and inhibition
    df_inh = df.loc[df['action'] == 'inhibition']

    # make separate source and target lists for activating and inhibiting
    sources_ut_a = [entry.split('.')[1] for entry in df_act['item_id_a']]
    targets_ut_a = [entry.split('.')[1] for entry in df_act['item_id_b']]
    sources_ut_i = [entry.split('.')[1] for entry in df_inh['item_id_a']]
    targets_ut_i = [entry.split('.')[1] for entry in df_inh['item_id_b']]

    # create edge with weight of 1
    edges_ut_a = zip(sources_ut_a, targets_ut_a, [1] * len(sources_ut_a))
    edges_ut_i = zip(sources_ut_i, targets_ut_i, [1] * len(sources_ut_i))
    edges_ut = edges_ut_a + edges_ut_i

    # create separate activating and inhibiting networks
    G_a = nx.MultiDiGraph()
    G_a.add_weighted_edges_from(edges_ut_a)

    G_i = nx.MultiDiGraph()
    G_i.add_weighted_edges_from(edges_ut_i)

    # add sign attribute
    edges_with_keys_a = G_a.edges(keys=True)
    signs_a = [1] * len(G_a.edges())
    edges_to_sign_a = dict(zip(edges_with_keys_a, signs_a))
    nx.set_edge_attributes(G_a, name='sign', values=edges_to_sign_a)

    edges_with_keys_i = G_i.edges(keys=True)
    signs_i = [-1] * len(G_i.edges())
    edges_to_sign_i = dict(zip(edges_with_keys_i, signs_i))
    nx.set_edge_attributes(G_i, name='sign', values=edges_to_sign_i)

    # combine two graphs
    G_str = nx.compose(G_a, G_i)

    # make quick translation list
    to_translate = list(set(zip(*edges_ut)[0] + zip(*edges_ut)[1]))

    # translate quick list and make translation dictionary
    mg = mygene.MyGeneInfo()
    mg_temp = mg.querymany(to_translate, scopes='ensemblprotein', fields='symbol')
    ensembl_list = [x['query'] for x in mg_temp]
    symbol_list = [x['symbol'] if 'symbol' in x.keys() else 'None' for x in mg_temp]
    ensembl_to_symbol = dict(zip(ensembl_list, symbol_list))

    # relabel nodes with symbols
    G_str = nx.relabel_nodes(G_str, ensembl_to_symbol)  # only keep the proteins that
    G_str.remove_node('None')

    return G_str


def filter_digraph(G,TF_list):

    if ((str(type(G)) == '<class \'networkx.classes.multidigraph.MultiDiGraph\'>') | (
                str(type(G)) == '<class \'networkx.classes.multigraph.MultiGraph\'>')):
        edges = G.out_edges(TF_list,keys = True, data = True)
        DG = nx.MultiDiGraph()
        DG.add_edges_from(edges)
    else:
        edges = G.out_edges(TF_list, data=True)
        DG = nx.DiGraph()
        DG.add_edges_from(edges)

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


