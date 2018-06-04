"""
-------------------------------------------
Author: Mikayla Webster (13webstermj@gmail.com)
        Brin Rosenthal (sbrosenthal@ucsd.edu)
Date: 10/13/17
-------------------------------------------
"""

import pandas as pd
import networkx as nx
import mygene # pip install mygene
import ndex2 # pip install ndex2
from ndex2.niceCXNetwork import NiceCXNetwork
import ndex2.client as nc


# ----------------------------- TRANSCRIPTION FACTOR -------------------------- #

def easy_load_TF_list(csv_filename = '../../TF_databases/TF_database_URA.csv', jaspar = True, TRED = True, ITFP = True, ENCODE = True, 
                      Neph2012 = True, TRRUST = True, Marbach2016 = True, species = 'human', gene_type = 'symbol'):
                      
    """
        Loads returns a list of transcription factors (TF) containing genes from the user specified databases.

        Args:
            csv_filename: String, filepath of where to find TF data ('../../TF_databases/TF_database_URA.csv' if you are using our 
                github directory structure)
                
            jaspar: Boolean, whether or not to include jaspar database in TF list
            TRED: Boolean, whether or not to include TRED database in TF list
            ITFP: Boolean, whether or not to include ITFP database in TF list
            ENCODE: Boolean, whether or not to include ENCODE database in TF list
            Neph2012: Boolean, whether or not to include Neph2012 database in TF list
            TRRUST: Boolean, whether or not to include TRRUST database in TF list
            Marbach2016: Boolean, whether or not to include Marbach2016 database in TF list
            
            species: String, specify either 'human' or 'mouse', should match the specs of your background database
            gene_type: String, specify either 'symbol' or 'entrez', should match the specs of your background database

        Returns: a list of transcription factors
    """
    
    df2 = pd.DataFrame.from_csv(csv_filename)
    TF_list = []
    
    if (jaspar == True):
        jaspar_list = list((df2.loc[df2['jaspar'] == 1])['jaspar'].index)
        TF_list.extend(jaspar_list)
        
    if (TRED == True):    
        TRED_list = list((df2.loc[df2['TRED'] == 1])['TRED'].index)
        TF_list.extend(TRED_list)
        
    if (ITFP == True):    
        ITFP_list = list((df2.loc[df2['ITFP'] == 1])['ITFP'].index)
        TF_list.extend(ITFP_list)
        
    if (ENCODE == True):    
        ENCODE_list = list((df2.loc[df2['ENCODE'] == 1])['ENCODE'].index)
        TF_list.extend(ENCODE_list)
        
    if (Neph2012 == True):   
        Neph_list = list((df2.loc[df2['Neph2012'] == 1])['Neph2012'].index)
        TF_list.extend(Neph_list)
        
    if (TRRUST == True):    
        TRRUST_list = list((df2.loc[df2['TRRUST'] == 1])['TRRUST'].index)
        TF_list.extend(TRRUST_list )
        
    if (Marbach2016 == True):    
        Marbach_list = list((df2.loc[df2['Marbach2016'] == 1])['Marbach2016'].index)
        TF_list.extend(Marbach_list)
        
    if (species == 'mouse'):
        TF_list = [x.capitalize() for x in TF_list]
    
    if (gene_type == 'entrez'):
        G_entrez = translate_gene_type(TF_list, 'symbol', 'entrezgene', species = species)
        TF_list = list(G_entrez.nodes())
        
    return list(set(TF_list))


# ------------------------- BACKGROUND NETWORK ------------------------------------ #


def load_small_STRING_to_digraph(filename, TF_list = []):

    """
        ** This is a depricated function that may not be compatible with the rest of our functions
    
        This function loads a small subset of the STRING database into a networkx digraph that can be used as input
        into the rest of our functions. This function can also filter the input database so that only the genes
        indicated by TF_list, and all of their out-going neighbors, will remain in the graph. Namely, the only
        source nodes left in the graph will be the genes in TF_list.

        ** Note that there is a loss of 110 our of about 40000 edges due to removal of multiedges **

        Args:
            filename: String, the path to the string file
            TF_list: List, the list of genes to filter the graph by

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
    DG_TF = keep_list_and_neighbors(DG_universe, TF_list)

    return DG_TF, DG_universe


# Use with STRING protein.actions files. This function is pretty exclusively for URA
def load_STRING_to_digraph(filename, TF_list = None, confidence_filter = 700, gene_type = 'symbol', species = 'human'):

    """
        This function loads a subset of the STRING database from either '9606.protein.actions.v10.5.txt' or
        '10090.protein.actions.v10.5.txt' into a networkx digraph that can be used as the background network for 
        the rest of our functions. IT ONLY KEEPS ACTIVATING AND INHIBITING EDGES FROM THE STRING NETWORK. This 
        function can also filter the input database so that only the genes indicated by TF_list, and all of their 
        out-going neighbors, will remain in the graph. Namely, the only source nodes left in the graph will be the 
        genes in TF_list. ** Intende for Upstream Regulator Analysis functions.

        Args:
            filename: String, the filepath to the STRING database file.
            TF_list: List, the list of genes to filter the graph by. If left as None, no filtering will happen.
            confidence_filter: A number between 0 and 1000, all interactions with confidece less than this number will be filtered out
            species: String, specify either 'human' or 'mouse', should match the specs of your differetial expression database
            gene_type: String, specify either 'symbol' or 'entrez', should match the specs of your differetial expression database

        Returns: 
            if TF_list is None, this funcitons returns:
                DG_universe: A Networkx DiGraph representation of the input STRING file
                
            if TF_list is supplied a list, this function returns:
                DG_TF:  A Networkx DiGraph representation of the input STRING file, filtered down to the sub-network of TF's and targets
                DG_universe: A Networkx DiGraph representation of the input STRING file

    """

    # read STRING file to dataframe
    df_full = pd.read_csv(filename, sep='\t')

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
	
    if TF_list == None:
        return DG_universe

    # filtered graph
    DG_TF = keep_list_and_neighbors(DG_universe, TF_list)

    return DG_TF, DG_universe
    
    
def load_ndex_from_server(UUID, relabel_node_field = None, filter_list = None):

    """
        This function loads an Ndex network as a networkx graph. This function can also filter the input database 
        so that only the genes indicated by filter_list, and all of their out-going neighbors, will remain in the 
        graph. This function is comparible to load_STRING_to_digraph().

        Args:
            UUID: String, the unique ID associated with the Ndex network you wish to load.
                  Example: 'e11c6684-5ac2-11e8-a4bf-0ac135e8bacf' <-- STRING_human_protein_actions_confidence_700_edges network
                           '76160faa-5d0f-11e8-a4bf-0ac135e8bacf' <-- STRING_mouse_protein_links_confidence_700_edges network
                           
            relabel_node_field: String, node attribute to relable the nodes in this graph with. Most Ndex graphs have a 'name' field
                  that contains the gene name you are looking for for each node. (Sometimes the nodes themselves are named weirdly).
                  In that case, set relabel_node_field = 'name'.
                  
            filter_list: List, the list of genes to filter the graph by. If None, no filtering will occur.

        Returns: 
            if filter_list is None, this funcitons returns:
                G: A Networkx representation of the input Ndex network
                
            if filter_list is supplied a list, this function returns:
                G_filtered:  A Networkx representation of the input Ndex network, filtered down to the sub-network of filter genes and filter gene targets
                G: A Networkx representation of the input Ndex network

    """
    
    # get nice object from server
    niceCx_from_server = ndex2.create_nice_cx_from_server(server='public.ndexbio.org', uuid=UUID)
    
    # convert to networkx graph
    G = niceCx_from_server.to_networkx()
    
    # rename the nodes to their gene symbol
    if relabel_node_field != None:
        node_id_to_name = nx.get_node_attributes(G, relabel_node_field)
        G = nx.relabel_nodes(G, node_id_to_name)
        #for (n,d) in G.nodes(data = True): del d[relabel_node_field] # delete the redundant field after
        
    # check to see if there is some weight field that needs to be converted to floats
    try:
        for n in G.edges(data = True):
            n[2]['weight'] = float(n[2]['weight'])
    except:
        pass

    # no filtering
    if filter_list == None:
        return G

    # filtered graph (say for TF's or DEG's)
    G_filtered = keep_list_and_neighbors(G, filter_list) # terribly named function. Actually does the opposite
    return G_filtered, G 

    
    
# use with STRING protein.links files  
def load_STRING_links(filename, confidence_filter = 700, species = 'human', translate_to = 'symbol'):

    """
        This function loads a subset of the STRING database from either '9606.protein.links.v10.5.txt' or
        '10090.protein.links.v10.5.txt' into a networkx graph that can be used as the background network for 
        the rest of our functions. **Intended for Heat Propogation and Clustering Analysis functions.

        Args:
            filename: String, the filepath to the STRING database file
            confidence_filter: A number between 0 and 1000, all interactions with confidece less than this number will be filtered out
            species: String, specify either 'human' or 'mouse', should match the specs of your differetial expression database
            translate_to: String, specify either 'symbol' or 'entrez', should match the specs of your differetial expression database

        Returns: 
            G_str: A Networkx representation of the input STRING file

    """
    
    # parse and load STRING protein links file
    string_df = pd.read_csv(filename, sep = ' ')
    string_df = string_df.loc[string_df['combined_score'] > confidence_filter]
                              
    # create the network
    G_str = nx.Graph()
    G_str.add_weighted_edges_from(zip(string_df['protein1'], string_df['protein2'], string_df['combined_score']))
    
    # extract list of ensemble gene names
    to_translate = [n[n.find('.') + 1:] for n in G_str.nodes()]
    G_str = nx.relabel_nodes(G_str, dict(zip(G_str.nodes(), to_translate)))

    # do the translation and relabeling using mygene
    G_str = translate_gene_type(to_translate, 'ensemblprotein', translate_to, G_str, species)
    return G_str
    
    
    

# --------------------- DEG LOAD FUNCTIONS ---------------------------#


def create_DEG_list(filename,
                    G = None,  # specify in order to add up-down info to graph
                    p_value_filter = 0.05,
                    p_value_or_adj = 'adj',  # filtering by p-value ('p') or adjusted p-value ('adj')
                    fold_change_filter = None,  # specify a number to filter by absolute (log) fold change
                    gene_type = 'symbol',  # 'symbol' or 'entrez'
                    gene_column_header = None,
                    p_value_column_header = None,
                    fold_change_column_header = None,
					sep = '\t',
                    return_full_values = False
                    ):

    """
        This function takes a differential expression data file and loads it into multiple dictionaries 
        that can be used by our other functions. Differential expression file must contain a gene name column, 
        a (log) fold change column, and an (adjusted) p-value column, each specified with a common column header.
        If your column header for any of these columns is not intuitive, specify the name of that column header.
        Each line in your input file should represent a different gene.

        Args:
            filename: String, filepath to the the differential expression input file
            G: Networkx graph, if specified will add (log) fold change data as node attribute 'updown' to G 
            p_value_filter: Float, typically a number between 0 and 1, the number to filter the adjusted p-vlaue by. Will remove all above this threshold
            p_value_or_adj: String, either 'p' or 'adj', that specifies whether we should use p-value or adjusted p-value information
            fold_change_filter: Float/Int, will filter out genes with fold change absolute value greater than this number
            gene_type: String, specify either 'symbol' or 'entrez', should match the specs of your background database
            gene_column_header: String, if your gene name column header is not intuitive, specify it here
            p_value_column_header: String, if your p-value column header is not intuitive, specify it here
            fold_change_column_header: String, if your fold change column header is not intuitive, specify it here
            sep: String, separating agent in your input file. Common exmaples include '\t' and ','
            return_full_values: Boolean, specifies whether to return only Differentially Expressed Genes' p-value and log fold
                information, or whether to return entire file's information.

        Returns:
            if parameter G is specified, return:
                DEG_list: list of remaining genes after all filtering is completed
                DG: Networkx Graph, a copy of G including 'updown' node attribute. Input graph G is not modified
                
            if parameter G is None, return:
                DEG_list: list of remaining genes after all filtering is completed
                DEG_to_pvalue: dictionary mapping gene symbol to (adjusted) p-value
                DEG_to_updown: dictionary mapping gene symbol to (log) fold change

    """

    df = pd.DataFrame.from_csv(filename, sep = sep)

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
    
    # save (log) fold change and p-value of entire graph, not just DEG's
    if return_full_values == True:
        gene_to_lfc = dict(zip(df[gene_column_header], df[fold_change_column_header]))
        gene_to_pvalue = dict(zip(df[gene_column_header], df[p_value_column_header]))

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
        
    # return early, because we don't want to filter down the pvalues and lfc 
    if return_full_values == True:
        return DEG_list, gene_to_pvalue, gene_to_lfc

    DEG_to_pvalue = dict(zip(DEG_list, df[p_value_column_header]))
    DEG_to_updown = dict(zip(DEG_list, df[fold_change_column_header]))

    # if a graph is provided, then add attribute to grpah rather than return as a dict
    if G != None:
        DG = add_node_attribute_from_dict(G, DEG_to_updown, attribute = 'updown')
        return DEG_list, DG
    
    return DEG_list, DEG_to_pvalue, DEG_to_updown



def create_DEG_full_graph(filename,
                    p_value_or_adj = 'adj',  # filtering by p-value ('p') or adjusted p-value ('adj')
                    gene_type = 'symbol',  # 'symbol' or 'entrez'
                    gene_column_header = None,
                    p_value_column_header = None,
                    fold_change_column_header = None,
					sep = '\t'
                    ):
                    
    """
        A special call to create_DEG_list() that loads the entire differetial expression file, and adds that information to
        a networkX graph.

        Args:
            filename: String, filepath to the the differential expression input file
            p_value_or_adj: String, either 'p' or 'adj', that specifies whether we should use p-value or adjusted p-value information
            gene_type: String, specify either 'symbol' or 'entrez', should match the specs of your background database
            gene_column_header: String, if your gene name column header is not intuitive, specify it here
            p_value_column_header: String, if your p-value column header is not intuitive, specify it here
            fold_change_column_header: String, if your fold change column header is not intuitive, specify it here
            sep: String, separating agent in your input file. Common exmaples include '\t' and ','

        Returns:
            DEG_full_graph: Networkx graph containing all genes in input file, and nodes attributes 'adj_p_value' and 'updown'
            DEG_to_pvalue: dictionary mapping gene symbol to (adjusted) p-value
            DEG_to_updown: dictionary mapping gene symbol to (log) fold change

    """

    DEG_full_list, DEG_to_pvalue, DEG_to_updown = create_DEG_list(filename,
                    p_value_filter = 1,
                    p_value_or_adj = p_value_or_adj,  # filtering by p-value ('p') or adjusted p-value ('adj')
                    gene_type = gene_type,  # 'symbol' or 'entrez'
                    gene_column_header = gene_column_header,
                    p_value_column_header = p_value_column_header,
                    fold_change_column_header = fold_change_column_header,
					sep = sep
                    )

    DEG_full_graph = nx.DiGraph()
    DEG_full_graph.add_nodes_from(DEG_full_list)
    DEG_full_graph = add_node_attribute_from_dict(DEG_full_graph, DEG_to_pvalue, attribute = 'adj_p_value')
    DEG_full_graph = add_node_attribute_from_dict(DEG_full_graph, DEG_to_updown, attribute = 'updown')

    return DEG_full_graph, DEG_to_pvalue, DEG_to_updown


# ------------------- HELPER FUNCTIONS ------------------------------ #


#def load_newline_sep_file(filename, column_header = 'genes'):
#
#    """
#        **Depricated. No longer in use.
#        Helper function to load TF files.
#    """
#
#    df = pd.read_csv(filename, names = [column_header], header = None)
#    return_list = list(df[column_header])
#    return return_list


def translate_gene_type(to_translate, before_gene_type, after_gene_type, G = None, species = 'human'):

    """
        Translates a list of genes using mygene. Returns those genes as a graph, where the nodes are the \
        translated genes.

        Args:
            to_translate: List, list of genes you would liek to translate.
            before_gene_type: String, gene naming convention of to_translate (symbol, entrez, etc.)
                Must match mygene naming conventions.
            after_gene_type: String, gene naming convention you wish to translate to. (symbol, entrez, etc...).
                Must match mygene naming conventions.
            G: Networkx Graph, if specified, G's nodes will be relabled using the new, translated genes. if left
                as None, a new networkx graph will be created and used.
            species: String, the species of the genes you are translating (ex. mouse, human, etc.) 
                Must match mygene naming conventions.

        Returns:
            G_after: A networkx graph whose nodes have been relabled with the translated genes. 

    """

    mg = mygene.MyGeneInfo()
    mg_temp = mg.querymany(to_translate, scopes = before_gene_type, fields = after_gene_type, species = species)
    before_list = [x['query'] for x in mg_temp]
    after_list = [str(x[after_gene_type]) if after_gene_type in x.keys() else 'None' for x in mg_temp]
    before_to_after = dict(zip(before_list, after_list))

    if G == None:
        G_before = nx.Graph()
        G_before.add_nodes_from(to_translate)
    else:
        G_before = G

    G_after = nx.relabel_nodes(G_before, before_to_after)  # only keep the proteins that
    if 'None' in G_after.nodes():
        G_after.remove_node('None')
    return G_after

    
def keep_list_and_neighbors(G, filter_list):

    """
        This is a helper function that removes all source nodes from graph G that are not in filter_list. The original
        graph is not modified.

        Args:
            G: Networkx graph, the graph to filter
            filter_list: List, the source nodes to keep

        Returns: A subset of the input networkx graph.

    """
    
    if type(G) == nx.classes.digraph.DiGraph:
        edges = G.out_edges(TF_list, data=True)
        DG = nx.DiGraph()
        DG.add_edges_from(edges)
        
    elif type(G) == nx.classes.multidigraph.MultiDiGraph:
        edges = G.out_edges(filter_list, keys = True, data = True)
        DG = nx.MultiDiGraph()
        DG.add_edges_from(edges)
        
    else:
        edges_to_keep = []
        [edges_to_keep.extend(nx.neighbors(G, n)) for n in filter_list]
        edges_to_keep.extend(filter_list)
        edges_to_keep = list(set(edges_to_keep))
        DG = G.subgraph(edges_to_keep)

    return DG



def add_node_attribute_from_dict(G, DEG_to_updown, attribute):

    """
        This function adds "updown" node attribute to nodes specified by DEG_to_updown. Does not modify 
        the input graph. Basically is a call to nx.set_node_attributes, but also handles some corner cases.

        Args:
            G: NetworkX graph, to add node attribute to.
            DEG_to_updown: Dict, maps gene symbol to up/down regulation information (output of create_DEG_list)

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

    """
        Helper function that determines if a header is in a DataFrame
    """

    try:
        df[header]
    except:
        return 0
    return 1


def try_message(df, header, message):

    """
        Helper function that returns an error message if a header is not in a DataFrame.
    """
    
    try:
        df[header]
    except:
        print(message)
        return -1
    return header


def check_gene_header(df, gene_column_header, gene_type):

    """
        Helper function that searches for common gene header names in a DataFrame. If none
        is found, will return an error saying so.
    """

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
            print('We could not find a gene column header in your file.')
            print('Please specify one with parameter \'gene_column_header\'.\n')
            return -1


def check_p_value_header(df, p_value_column_header, p_value_or_adj):

    """
        Helper function that searches for common p-value or adjusted p-value header names in a DataFrame. If none
        is found, will return an error saying so.
    """

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
            print('\"p_value_or_adj\" argument can be only \'adj\' or \'p\'. It is currently set to \'' + p_value_or_adj + '\'.\n')

        for header in common_p_value_headers:
            if try_or(df, header) == 1:  # if we find a matching header
                return header

        if p_value_column_header == None:  # if we did not find a matching header
            print('We could not find a p-value column header in your file.')
            print('Please specify one with parameter \'p_value_column_header\'.\n')
            return -1



def check_fold_change_header(df, fold_change_column_header):

    """
        Helper function that searches for common fold change or log fold change header names in a DataFrame. If none
        is found, will return an error saying so.
    """

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
            print('We could not find a fold change column header in your file.')
            print('Please specify one with parameter \'fold_change_column_header\'.\n')
            return -1


