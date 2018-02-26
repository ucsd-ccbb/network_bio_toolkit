"""
-------------------------------------------
Author: Mikayla Webster (m1webste@ucsd.edu)
Date: 10/13/17
-------------------------------------------
"""

import pandas as pd
import networkx as nx
import mygene


# ----------------------------- TRANSCRIPTION FACTOR -------------------------- #


def load_slowkow(filename_list=['../../TF_databases/slowkow_databases/TRED_TF.txt',
                                   '../../TF_databases/slowkow_databases/ITFP_TF.txt',
                                   '../../TF_databases/slowkow_databases/ENCODE_TF.txt',
                                   '../../TF_databases/slowkow_databases/Neph2012_TF.txt',
                                   '../../TF_databases/slowkow_databases/TRRUST_TF.txt',
                                   '../../TF_databases/slowkow_databases/Marbach2016_TF.txt']):

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
    for filename in filename_list:
        to_add = load_newline_sep_file(filename)
        return_list.extend([x.upper() for x in to_add])

    # remove duplicates
    return set(return_list)




def load_jaspar(filename = 'jaspar_genereg_matrix.txt'):

    """
        This function loads The 2016 version of Jaspar's TF database, which can be dowloaded as a .txt file from
        http://jaspar2016.genereg.net/. At the bottom of the page, navigate to the "Download" button. Open the
        "database" file and download "MATRIX.txt".
        Args:
            filename: string, the path to the jaspar "MATRIX.txt" file
        Returns: A list of all unique transcription factors from the file provided
    """

    # parse jaspar file
    jasp_df = pd.read_csv(filename, sep="\t", header=None, names=['col1', 'col2', 'col3', 'col4', 'tf_genes'])

    # return transcription factors with ALL CAPS names
    return list(set(jasp_df['tf_genes'].str.upper()))



def easy_load_TF_list(slowkow_bool=True,
                    slowkow_files=['../../TF_databases/slowkow_databases/TRED_TF.txt',
                                   '../../TF_databases/slowkow_databases/ITFP_TF.txt',
                                   '../../TF_databases/slowkow_databases/ENCODE_TF.txt',
                                   '../../TF_databases/slowkow_databases/Neph2012_TF.txt',
                                   '../../TF_databases/slowkow_databases/TRRUST_TF.txt',
                                   '../../TF_databases/slowkow_databases/Marbach2016_TF.txt'],
                    jaspar_bool=True,
                    jaspar_file="../../TF_databases/jaspar_genereg_matrix.txt",
                    gene_type = "symbol",
                    species = 'human'):
    """
        This function loads The 2016 version of Jaspar's TF database, which can be dowloaded as a .txt file from
        http://jaspar2016.genereg.net/. At the bottom of the page, navigate to the "Download" button. Open the
        "database" file and download "MATRIX.txt".

        Args:
            slowkow_bool: boolean, whether or not to include the slowkow databases
            slowkow_files: list or strings, a list of input file path names to the load_slowkow method
            slowkow_sep: string, the type of deliminator used in the slowkow_files to separate the transcription factors
            jaspar_bool: boolean, whether or not to include the jaspar database
            jaspar_file: string, the file path where to find the jaspar database

        Returns: A list of all unique transcription factors from the files provided

    """

    TF_list = []

    if slowkow_bool == True:
        slowkow_TFs = load_slowkow(slowkow_files)
        TF_list.extend(slowkow_TFs)

    if jaspar_bool == True:
        jaspar_TFs = load_jaspar(jaspar_file)
        TF_list.extend(jaspar_TFs)

    if gene_type == 'entrez':
        G_entrez = translate_gene_type(TF_list, 'symbol', 'entrezgene', species = species)
        TF_list = list(G_entrez.nodes())

    return list(set(TF_list))


# ------------------------- BACKGROUND NETWORK ------------------------------------ #


def load_small_STRING_to_digraph(filename, TF_list = []):

    """
        This function loads a small subset of the STRING database into a networkx digraph that can be used as input
        into the rest of our functions. This function can also filter the input database so that only the genes
        indicated by TF_list, and all of their out-going neighbors, will remain in the graph. Namely, the only
        source nodes left in the graph will be the genes in TF_list.

        ** Note that there is a loss of 110 our of about 40000 edges due to removal of multiedges **

        Args:
            filename: string, the path to the string file
            TF_list: the list of genes to filter the graph by

        Returns: DG_universe: A Networkx DiGraph representation, using all-caps gene symbol, of the input STRING file
				 DG_TF:  A Networkx DiGraph representation of the input STRING file, filtered down to the sub-network of TF's and targets

    """

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
    DG_universe = nx.from_pandas_dataframe(STRING_DF, 'source', 'target', ['sign', 'weight'], create_using=nx.DiGraph())

    # filtered graph
    DG_TF = remove_source_nodes_from_list(DG_universe, TF_list)

    return DG_TF, DG_universe



def load_STRING_to_digraph(filename, TF_list, confidence_filter=400, gene_type = "symbol", species = 'human'):

    """
        This function loads a subset of the STRING database from the file "9606.protein.actions.v10.5.txt"
        into a networkx digraph that can be used as input into the rest of our functions. This function can
        also filter the input database so that only the genes indicated by TF_list, and all of their out-going
        neighbors, will remain in the graph. Namely, the only source nodes left in the graph will be the genes in TF_list.

        Args:
            filename: string, the path to the string file
            confidence_filter: A number between 0 and 1000, all interactions with confidece less than this number
                will be filtered out
            TF_list: the list of genes to filter the graph by

        Returns: DG_universe: A Networkx DiGraph representation, using all-caps gene symbol, of the input STRING file
				 DG_TF:  A Networkx DiGraph representation of the input STRING file, filtered down to the sub-network of TF's and targets

    """

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

    # decide which "language" to translate_gene_type to
    before_gene_type = 'ensemblprotein'
    if gene_type == 'entrez':
        after_gene_type = 'entrezgene'
    else:
        after_gene_type = 'symbol'

    # do the translation
    DG_universe = translate_gene_type(to_translate, before_gene_type, after_gene_type, G_str, species)

    # filtered graph
    DG_TF = remove_source_nodes_from_list(DG_universe, TF_list)

    return DG_TF, DG_universe


# --------------------- DEG LOAD FUNCTIONS ---------------------------#


def create_DEG_list(filename,
                    G=None,  # specify in order to add up-down info to graph

                    p_value_filter=0.05,
                    p_value_or_adj='adj',  # filtering by p-value ('p') or adjusted p-value ('adj')

                    fold_change_filter=None,  # specify a number to filter by absolute (log) fold change
                    gene_type='symbol',  # 'symbol' or 'entrez'

                    gene_column_header=None,
                    p_value_column_header=None,
                    fold_change_column_header=None
                    ):

    """
        This function takes a standard input file representation of a list of differentially expressed genes and
        loads it into multiple dictionaries that can be used by our functions. Our standard input file must be
        a tab separated list, where each row represents a gene. This file must contain column headers "adj_p_value"
        (adjusted p-value), "gene_symbol", and "fold_change" (the fold change or log fold change).

        Args:
            filename: the standard input file
            p_value_filter: typically a number between 0 and 1, the number to filter the adjusted p-vlaue by. Will remove all above this threshold
			fold_change_filter: Will remove all below this threshold

        Returns:
            DEG_list: list of gene symbols as strings
            DEG_to_pvalue: dictionary mapping DEG gene symbol to adjusted p-value
            DEG_to_updown: dictionary mapping DEG gene symbol to (log) fold change

    """

    df = pd.DataFrame.from_csv(filename, sep='\t')

    # check to make sure we know which column headers to use
    gene_column_header = check_gene_header(df, gene_column_header, gene_type)
    p_value_column_header = check_p_value_header(df, p_value_column_header, p_value_or_adj)
    fold_change_column_header = check_fold_change_header(df, fold_change_column_header)

    if (gene_column_header == -1) | (p_value_column_header == -1) | (fold_change_column_header == -1):
        if G == None:
            return None, None
        return None, None, None

    # remove duplicate lines for same gene symbol, just use first occurance
    df.drop_duplicates(subset=[gene_column_header], keep='first', inplace=True)

    # filter by p-value cut-off
    df = df.loc[df[p_value_column_header] < p_value_filter]

    # filter by (log) fold change cut off if applicable
    if fold_change_filter != None:
        df = df.loc[abs(df[fold_change_column_header]) > fold_change_filter]

    DEG_list = df[gene_column_header]

    DEG_list = list(DEG_list)
    if (type(DEG_list[0]) != str):
        DEG_list_temp = [str(int(x)) for x in DEG_list if (str(x) != 'nan')]
        DEG_list = DEG_list_temp

    DEG_to_pvalue = dict(zip(DEG_list, df[p_value_column_header]))
    DEG_to_updown = dict(zip(DEG_list, df[fold_change_column_header]))

    if G != None:
        DG = add_node_attribute_from_dict(G, DEG_to_updown, attribute='updown')
        return DEG_list, DG

    return DEG_list, DEG_to_pvalue, DEG_to_updown



def create_DEG_full_graph(filename,
                    p_value_or_adj='adj',  # filtering by p-value ('p') or adjusted p-value ('adj')
                    gene_type='symbol',  # 'symbol' or 'entrez'
                    gene_column_header=None,
                    p_value_column_header=None,
                    fold_change_column_header=None
                    ):

    DEG_full_list, DEG_to_pvalue, DEG_to_updown = create_DEG_list(filename,
                    p_value_filter=1,
                    p_value_or_adj=p_value_or_adj,  # filtering by p-value ('p') or adjusted p-value ('adj')
                    gene_type=gene_type,  # 'symbol' or 'entrez'
                    gene_column_header=gene_column_header,
                    p_value_column_header=p_value_column_header,
                    fold_change_column_header=fold_change_column_header
                    )

    DEG_full_graph = nx.DiGraph()
    DEG_full_graph.add_nodes_from(DEG_full_list)
    DEG_full_graph = add_node_attribute_from_dict(DEG_full_graph, DEG_to_pvalue, attribute='adj_p_value')
    DEG_full_graph = add_node_attribute_from_dict(DEG_full_graph, DEG_to_updown, attribute='updown')


    return DEG_full_graph, DEG_to_pvalue, DEG_to_updown


# ------------------- HELPER FUNCTIONS ------------------------------ #


def load_newline_sep_file(filename, column_header = 'genes'):
    df = pd.read_csv(filename, names=[column_header], header=None)
    return_list = list(df[column_header])
    return return_list


def translate_gene_type(to_translate, before_gene_type, after_gene_type, G = None, species = 'human'):

    mg = mygene.MyGeneInfo()
    mg_temp = mg.querymany(to_translate, scopes=before_gene_type, fields=after_gene_type, species = species)
    before_list = [x['query'] for x in mg_temp]
    after_list = [str(x[after_gene_type]) if after_gene_type in x.keys() else 'None' for x in mg_temp]
    before_to_after = dict(zip(before_list, after_list))

    if G == None:
        G_before = nx.Graph()
        G_before.add_nodes_from(to_translate)
    else:
        G_before = G

    G_after = nx.relabel_nodes(G_before, before_to_after)  # only keep the proteins that
    G_after.remove_node('None')
    return G_after


def remove_source_nodes_from_list(G, TF_list):

    """
        This is a helper function that removes all source nodes from graph G that are not in TF_list. The original
        graph is not modified.

        Args:
            G: The NetworkX graph to filter
            TF_list: the source nodes to keep

        Returns: A Networkx graph representation of the STRING database file

    """

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



def add_node_attribute_from_dict(G, DEG_to_updown, attribute):

    """
        This function adds "updown" node attribute to nodes specified by DEG_to_updown

        Args:
            G: DiGraph, a networkX directed graph
            DEG_to_updown: a dictionary that maps gene symbol to up/down regulation information (output of create_DEG_list)

        Returns: A networkX graph, a copy of the input graph with added node attribute "updown" to applicable genes
    """

    # don't want to add updowns to original graph
    if str(type(G)) == '<class \'networkx.classes.graph.Graph\'>':
        DG = nx.Graph(G)
    elif str(type(G)) == '<class \'networkx.classes.digraph.DiGraph\'>':
        DG = nx.DiGraph(G)
    elif str(type(G)) == '<class \'networkx.classes.multigraph.MultiGraph\'>':
        DG = nx.MultiGraph(G)
    elif str(type(G)) == '<class \'networkx.classes.multidigraph.MultiDiGraph\'>':
        DG = nx.MultiDiGraph(G)
    else:
        return -1

    # get all the differentially expressed genes in DG
    DEG_in_DG = set(DG.nodes()) & set(DEG_to_updown.keys())

    # add node attribute to each node in DG if it exists, otherwise set to zero
    zero_dict = dict(zip(DG.nodes(), [0] * len(DG.nodes())))
    for gene in DEG_in_DG:
        zero_dict[gene] = DEG_to_updown[gene]
    nx.set_node_attributes(DG, attribute, zero_dict)

    return DG


def try_or(df, header):
    try:
        df[header]
    except:
        return 0
    return 1


def try_message(df, header, message):
    try:
        df[header]
    except:
        print message
        return -1
    return header


def check_gene_header(df, gene_column_header, gene_type):

    # if the user specified a column header
    if gene_column_header != None:
        error_message = "The gene column header you specified does not exist in your file.\n"
        gene_column_header = try_message(df, gene_column_header, error_message)
        return gene_column_header

    # or check if it is going by a stardardly excepted format
    else:
        if gene_type == 'symbol':
            common_gene_headers = ['gene', 'GENE', 'symbol', 'SYMBOL', 'genesymbol', 'GENESYMBOL',
                                   'gene_symbol', 'GENE_SYMBOL', 'gene symbol', 'GENE SYMBOL']
        elif gene_type == 'entrez':
            common_gene_headers = ['gene', 'GENE', 'entrez', 'ENTREZ', 'entrezid', 'ENTREZID',
                                   'entrez_id', 'ENTREZ_ID', 'entrez id', 'ENTREZ ID']

        for header in common_gene_headers:
            if try_or(df, header) == 1:  # if we find a matching header
                return header

        if gene_column_header == None:  # if we did not find a matching header
            print 'We could not find a gene column header in your file.'
            print 'Please specify one with parameter \'gene_column_header\'.\n'
            return -1


def check_p_value_header(df, p_value_column_header, p_value_or_adj):

    # if the user specified a column header
    if p_value_column_header != None:
        error_message = "The p-value column header you specified does not exist in your file.\n"
        p_value_column_header = try_message(df, p_value_column_header, error_message)
        return p_value_column_header

    # or check if it is going by a stardardly excepted format
    else:
        if p_value_or_adj == 'adj':
            common_p_value_headers = ['adj.P.Val', 'padj', 'adj_p_value', 'ADJ_P_VALUE',
                                      'adjP', 'ADJP', 'adjp', 'adj_p_val', 'ADJ_P_VAL', 'adjusted p-value']
        elif p_value_or_adj == 'p':
            common_p_value_headers = ['P.Value', 'pvalue', 'PVALUE', 'p-value', 'P-VALUE',
                                      'p.value', 'P.VALUE', 'pval', 'PVAL', 'p_val', 'P_VAL']
        else:
            print '\"p_value_or_adj\" argument can be only \'adj\' or \'p\'. It is currently set to \'' + p_value_or_adj + '\'.\n'

        for header in common_p_value_headers:
            if try_or(df, header) == 1:  # if we find a matching header
                return header

        if p_value_column_header == None:  # if we did not find a matching header
            print 'We could not find a p-value column header in your file.'
            print 'Please specify one with parameter \'p_value_column_header\'.\n'
            return -1



def check_fold_change_header(df, fold_change_column_header):

    # if the user specified a column header
    if fold_change_column_header != None:
        error_message = "The fold change column header you specified does not exist in your file.\n"
        fold_change_column_header = try_message(df, fold_change_column_header, error_message)
        return fold_change_column_header

    # or check if it is going by a stardardly excepted format
    else:
        common_fold_change_headers = ['logFC', 'log2FoldChange', 'log_fold_change', 'LOG_FOLD_CHANGE',
                                      'log-fold-change', 'LOG-FOLD-CHANGE', 'lfc', 'LFC', 'log fold change'
                                                                                          'LOG FOLD CHANGE', 'fldchg',
                                      'FLDCHG', 'fc', 'FC', 'fold_change',
                                      'FOLD_CHANGE', 'fold change', 'FOLD CHANGE', 'foldchange', 'FOLDCHANGE',
                                      'fold-change', 'FOLD-CHANGE']

        for header in common_fold_change_headers:
            if try_or(df, header) == 1:  # if we find a matching header
                return header

        if fold_change_column_header == None:  # if we did not find a matching header
            print 'We could not find a fold change column header in your file.'
            print 'Please specify one with parameter \'fold_change_column_header\'.\n'
            return -1


