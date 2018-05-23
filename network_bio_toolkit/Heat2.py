"""
-------------------------------------------
Author: Mikayla Webster (13webstermj@gmail.com)
Date: 4/6/18
-------------------------------------------
"""

import create_graph
import heat_and_cluster
reload(heat_and_cluster)
import visJS2jupyter.visualizations as visualizations # pip install visJS2jupyter

#for local testing
#import sys
#sys.path.append('../../../visJS2jupyter/visJS2jupyter')
#import visualizations
#import visJS_module

import networkx as nx
import pandas as pd

class Heat:

# --------------- INITIATION ------------------------------------- #

    def __init__(self, gene_type = 'symbol', species = 'human'):
        
        # User has to provide a gene type
        if (gene_type != 'entrez') & (gene_type != 'symbol'):
            raise ValueError('Please specify either \"entrez\" or \"symbol\" as your initialization argument.\n')

        # User has to provide a gene type
        if (species != 'human') & (species != 'mouse'):
            raise ValueError('Please specify either \"entrez\" or \"symbol\" as your initialization argument.\n')

        # instantiate all of our instance variables
        self.gene_type = gene_type
        self.species = species
        self.DEG_list = None
        self.return_entire_lf = None
        self.DG_universe = None
        self.Wprime = None
        
        # map string to actual instance variable
        self.string_to_item = {}
        self.string_to_item['gene_type'] = self.gene_type
        self.string_to_item['species'] = self.species
        self.string_to_item['DEG_list'] = self.DEG_list
        self.string_to_item['return_entire_lf'] = self.return_entire_lf
        self.string_to_item['DG_universe'] = self.DG_universe
        self.string_to_item['Wprime'] = self.Wprime
        
        # instanciating error message dict
        self.item_to_message = {}
        self.item_to_message['gene_type'] = 'No gene type currently on file. Please create a new instance of Heat using '\
                     + 'Heat(gene_type, species)'
        self.item_to_message['species'] = 'No species declaration currently on file. Please create a new instance of Heat using ' \
                     + 'Heat(gene_type, species)'
        self.item_to_message['DEG_list'] = 'No differentially expressed gene list currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.create_DEG_list\n' \
                     + 'Or assign your own using Heat_instance.DEG_list\n'
        self.item_to_message['return_entire_lf'] = 'No fold change information currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.create_DEG_list\n' \
                     + 'Or assign your own using Heat_instance.return_entire_lf\n'
        self.item_to_message['DG_universe'] = 'No background network currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.load_STRING_to_digraph\n' \
                     + 'Or assign your own using Heat_instance.DG_universe\n'
        self.item_to_message['Wprime'] = 'No adjacency matrix currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.normalized_adj_matrix()\n' \
                     + 'Or assign your own using Heat_instance.Wprime\n'

            
# ----------------------------- GETTERS AND SETTERS -------------------------- #
# So our users don't accidentally modify an instance variable without realizing

    # to ask to see a variable
    def get(self, item):

        # check that the argument is valid, and update our dictionary
        exists = self.check_exists(item)
        if exists == False:
            return None

        # check if the user input is valid
        to_return = self.string_to_item[item]

        # check if the type of that variable is a "primitive" type
        return_deep_copy = True
        basic_python_types = [list, dict, tuple, float, int, str, bool]
        for py_type in basic_python_types:
            if type(to_return) == py_type:
                return_deep_copy = False
                break

        # if it is a primitive type, just return it. Else return a deep copy
        if return_deep_copy == True:
            return copy.deepcopy(to_return)
        else:
            return to_return


    # to set the value of a variable
    def set(self, item, value):

        self.check_exists(item)

        # check if the type of that variable is a "primitive" type
        return_deep_copy = True
        basic_python_types = [list, dict, tuple, float, int, str, bool]
        for py_type in basic_python_types:
            if type(value) == py_type:
                return_deep_copy = False
                break

        # if it is not a "primitive" type, copy it, just in case
        if return_deep_copy == True:
            value = copy.deepcopy(value)

        # map string to actual instance variable, then set that instance variable
        if item == 'gene_type': self.gene_type = value
        elif item == 'species': self.species = value
        elif item == 'DEG_list': self.DEG_list = value
        elif item == 'return_entire_lf': self.return_entire_lf = value
        elif item == 'DG_universe': self.DG_universe = value
        elif item == 'Wprime': self.Wprime = value
        
        else:
            print ('The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
            + '- gene_type\n' \
            + '- species\n' \
            + '- DEG_list\n' \
            + '- return_entire_lf\n' \
            + '- DG_universe\n' \
            + '- Wprime\n\n')


#----------------------- LOAD NETWORK FUNCTIONS ---------------------------------#
            
    def create_DEG_list(self, filename,
                        p_value_filter=0.05,
                        p_value_or_adj='adj',  # filtering by p-value ('p') or adjusted p-value ('adj')
                        fold_change_filter=None,  # specify a number to filter by absolute (log) fold change
                        gene_column_header=None,
                        p_value_column_header=None,
                        fold_change_column_header=None,
                        sep = '\t'):

        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        # create the DEG list with specified cut-offs
        DEG_list, DEG_to_pvalue, return_entire_lf = create_graph.create_DEG_list(filename, None, p_value_filter, p_value_or_adj,
                fold_change_filter, self.gene_type, gene_column_header, p_value_column_header, fold_change_column_header, sep,
                return_entire_lf = True)
        self.DEG_list = DEG_list
        self.return_entire_lf = return_entire_lf
        
        
    def load_STRING_to_digraph(self, filename, confidence_filter=400):

        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        DG_universe = create_graph.load_STRING_to_digraph(filename, None, confidence_filter, self.gene_type, self.species)
        self.DG_universe = DG_universe
        
    def load_ndex_from_server(self, UUID, relabel_node_field = None):
    
        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return
    
        DG_universe = create_graph.load_ndex_from_server(UUID, relabel_node_field, None)
        self.DG_universe = DG_universe
        
    def load_STRING_links(self, filename, confidence_filter = 700):
    
        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return
    
        # load STRING
        self.DG_universe = create_graph.load_STRING_links(filename, confidence_filter, self.species, self.gene_type)
        
        
#------------------------- Heat Propagation --------------------------------#

    def normalized_adj_matrix(self):
    
        # make sure user has run all prerequisites
        if self.check_exists('DG_universe') == False: return None
        self.Wprime = visualizations.normalized_adj_matrix(nx.Graph(self.DG_universe))
        
   
    def draw_heat_prop(self, 
                        num_nodes = 200,
                        edge_width = 2,
                        node_size_multiplier = 5,
                        largest_connected_component = True,
                        physics_enabled = True,
                        node_font_size = 40,
                        **kwargs):
    
        # make sure user has run all prerequisites
        for item in ['DG_universe', 'DEG_list', 'Wprime']:
            if self.check_exists(item) == False:
                return

        G_heat = nx.Graph(self.DG_universe)
        seed_nodes = [n for n in self.DEG_list if n in self.DG_universe]
        
        return visualizations.draw_heat_prop(G_heat, seed_nodes,
                                            Wprime = self.Wprime,
                                            num_nodes = num_nodes,
                                            highlight_nodes = seed_nodes,
                                            edge_width = edge_width,
                                            node_size_multiplier = node_size_multiplier,
                                            largest_connected_component = largest_connected_component,
                                            physics_enabled = physics_enabled,
                                            node_font_size = node_font_size,
                                            **kwargs)
                                            
                                            
#-------------------------- CLUSTERING -----------------------------------#

    def draw_clustering(self,
                    rad_positions = True,
                    k = None,
                    largest_connected_component = False,
                    alpha = 0.5,
                    num_its = 20,
                    num_top_genes = 200,
                    cluster_size_cut_off = 0,
                    remove_stray_nodes = False,
                    r = 2.0,
                    x_offset = 2,
                    y_offset = 2,
                    node_spacing = 500,
                    node_size_multiplier = 10,
                    physics_enabled = False,
                    node_font_size = 40,
                    graph_id = 3,
                    **kwargs
               ):

        # make sure user has run all prerequisites
        for item in ['DEG_list', 'DG_universe', 'Wprime']:
            if self.check_exists(item) == False:
                return
                
        seed_nodes = [n for n in self.DEG_list if n in self.DG_universe]

        return heat_and_cluster.draw_clustering(self.DG_universe, seed_nodes,
                    rad_positions = rad_positions,
                    Wprime = self.Wprime,
                    k = k,
                    largest_connected_component = largest_connected_component,
                    alpha = alpha,
                    num_its = num_its,
                    num_top_genes = num_top_genes,
                    cluster_size_cut_off = cluster_size_cut_off,
                    remove_stray_nodes = remove_stray_nodes,
                    r = r,
                    x_offset = x_offset,
                    y_offset = y_offset,
                    node_spacing = node_spacing,
                    node_size_multiplier = node_size_multiplier,
                    physics_enabled = physics_enabled,
                    node_font_size = node_font_size,
                    graph_id = graph_id,
                    return_entire_lf = self.return_entire_lf,
                    **kwargs
                    )
                    
                    


#------------------ HELPER FUNCTIONS ------------------------------------#  

 
    # item must be string version
    def check_exists(self, item):

        # re-map it so it stays up to date
        self.string_to_item = {}
        self.string_to_item['gene_type'] = self.gene_type
        self.string_to_item['species'] = self.species
        self.string_to_item['DEG_list'] = self.DEG_list
        self.string_to_item['return_entire_lf'] = self.return_entire_lf
        self.string_to_item['DG_universe'] = self.DG_universe
        self.string_to_item['Wprime'] = self.Wprime

        try:
            if (type(self.string_to_item[item]) == type(None)):
                print(self.item_to_message[item])
                return False
        except:
            print('The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
            + '- gene_type\n' \
            + '- species\n' \
            + '- DEG_list\n' \
            + '- return_entire_lf\n' \
            + '- DG_universe\n' \
            + '- Wprime\n\n')
            return False
        return True