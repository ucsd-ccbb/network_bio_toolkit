"""
-------------------------------------------
Author: Mikayla Webster (13webstermj@gmail.com)
Date: 4/6/18
-------------------------------------------
"""
# our modules
import create_graph
import stat_analysis
import heat_and_cluster
import visJS2jupyter.visualizations as visualizations # pip install visJS2jupyter

# common packages, most liekly already installed
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy

# reloading, for testing
reload(create_graph)
reload(stat_analysis)
reload(heat_and_cluster)

#for local testing
#import sys
#sys.path.append('../../../visJS2jupyter/visJS2jupyter')
#import visualizations
#import visJS_module

class Heat:

# --------------- INITIATION ------------------------------------- #

    def __init__(self, gene_type = 'symbol', species = 'human'):
    
        """
            Initiates Heat class instance, initiating all instance variables to None.
        """
        
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
        self.node_to_lfc = None
        self.node_to_pvalue = None
        self.DG_universe = None
        self.Wprime = None
        self.node_to_cluster = None
        
        # map string to actual instance variable
        self.string_to_item = {}
        self.string_to_item['gene_type'] = self.gene_type
        self.string_to_item['species'] = self.species
        self.string_to_item['DEG_list'] = self.DEG_list
        self.string_to_item['node_to_lfc'] = self.node_to_lfc
        self.string_to_item['node_to_pvalue'] = self.node_to_pvalue
        self.string_to_item['DG_universe'] = self.DG_universe
        self.string_to_item['Wprime'] = self.Wprime
        self.string_to_item['node_to_cluster'] = self.node_to_cluster
        
        # instanciating error message dict
        self.item_to_message = {}
        self.item_to_message['gene_type'] = 'No gene type currently on file. Please create a new instance of Heat using '\
                     + 'Heat(gene_type, species)'
        self.item_to_message['species'] = 'No species declaration currently on file. Please create a new instance of Heat using ' \
                     + 'Heat(gene_type, species)'
        self.item_to_message['DEG_list'] = 'No differentially expressed gene list currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.create_DEG_list\n' \
                     + 'Or assign your own using Heat_instance.DEG_list\n'
        self.item_to_message['node_to_lfc'] = 'No fold change information currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.create_DEG_list\n' \
                     + 'Or assign your own using Heat_instance.node_to_lfc\n'
        self.item_to_message['node_to_pvalue'] = 'No p-value information currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.create_DEG_list\n' \
                     + 'Or assign your own using Heat_instance.pvalue\n'
        self.item_to_message['DG_universe'] = 'No background network currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.load_STRING_to_digraph\n' \
                     + 'Or assign your own using Heat_instance.DG_universe\n'
        self.item_to_message['Wprime'] = 'No adjacency matrix currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.normalized_adj_matrix()\n' \
                     + 'Or assign your own using Heat_instance.Wprime\n'
        self.item_to_message['node_to_cluster'] = 'No cluster dictionary currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.draw_clustering()\n' \
                     + 'Or assign your own using Heat_instance.node_to_cluster\n'
                     
                     
    def copyHeat(self):
    
        """
            Creates a deep copy of this analysis by initiating a new Upstream instance
            and (deep) copying all instance variables to this new instance.
        """

        to_return = Heat(self.gene_type, self.species)
        to_return.gene_type = self.gene_type
        to_return.species = self.species
        to_return.DEG_list = self.DEG_list
        to_return.node_to_lfc = self.node_to_lfc
        to_return.node_to_pvalue = self.node_to_pvalue
        to_return.DG_universe = copy.deepcopy(self.DG_universe)
        to_return.Wprime = copy.deepCopy(to_return.Wprime)
        to_return.node_to_cluster = self.node_to_cluster
        to_return.string_to_item = self.string_to_item
        to_return.item_to_message = self.item_to_message

        return to_return
        
            
    def __repr__(self):
        return self.create_to_string()
        

    def __str__(self):
        return self.create_to_string()
        

    def create_to_string(self):
    
        """
            Creates a string representation of this object by printing whether each of its
            instance variables has been instanciated or not. If an instance variable has been 
            instanciated, it will print the type of that variable.
        """
    
        for item in ['gene_type', 'species', 'DEG_list', 'node_to_lfc', 'node_to_pvalue', 
                     'DG_universe', 'Wprime', 'node_to_cluster']:
            print (item + ': ')
            exists = self.check_exists(item)
            if exists == True:
                print str(type(self.string_to_item[item])) + '\n'
        return '\n'

            
# ----------------------------- GETTERS AND SETTERS -------------------------- #

    # to ask to see a variable
    def get(self, item):
    
        """
            Use this function to access instance variable values, but not the actual instance variable
            references. (aka, if you modify when this function returns, it will not modify the stored
            version of that variable.)
            
            Arg:
                item: String, name of the variable you'd like to access. (Look at heat_instance.item_to_message 
                      if you are unsure.
        """

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
    
        """
            Use this function to update the value of an instance variable.
            
            Arg:
                item: String, name of the variable you'd like to modify. (Look at heat_instance.item_to_message 
                      if you are unsure.
        """

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
        elif item == 'node_to_lfc': self.node_to_lfc = value
        elif item == 'node_to_pvalue': self.node_to_pvalue = value
        elif item == 'DG_universe': self.DG_universe = value
        elif item == 'Wprime': self.Wprime = value
        elif item == 'node_to_cluster': self.node_to_cluster = value
        
        else:
            print ('The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
            + '- gene_type\n' \
            + '- species\n' \
            + '- DEG_list\n' \
            + '- node_to_lfc\n' \
            + '- node_to_pvalue\n' \
            + '- DG_universe\n' \
            + '- Wprime\n' \
            + '- node_to_cluster\n\n')    
            


#----------------------- LOAD NETWORK FUNCTIONS ---------------------------------#
            
    def create_DEG_list(self, filename,
                        p_value_filter=0.05,
                        p_value_or_adj='adj',  # filtering by p-value ('p') or adjusted p-value ('adj')
                        fold_change_filter=None,  # specify a number to filter by absolute (log) fold change
                        gene_column_header=None,
                        p_value_column_header=None,
                        fold_change_column_header=None,
                        sep = '\t'):
                        
        """
            Use this function to access instance variable values, but not the actual instance variable
            references. (aka, if you modify when this function returns, it will not modify the stored
            version of that variable.)
            
            Args:
                item: String, name of the variable you'd like to access. (Look at heat_instance.item_to_message 
                      if you are unsure.
        """

        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        # create the DEG list with specified cut-offs
        DEG_list, node_to_pvalue, node_to_lfc = create_graph.create_DEG_list(filename, None, p_value_filter, p_value_or_adj,
                fold_change_filter, self.gene_type, gene_column_header, p_value_column_header, fold_change_column_header, sep,
                return_full_values = True)
        self.DEG_list = DEG_list
        self.node_to_lfc = node_to_lfc
        self.node_to_pvalue = node_to_pvalue
        
        
    def load_STRING_to_digraph(self, filename, confidence_filter = 400):
    
        """
            This function loads a subset of the STRING database from either '9606.protein.actions.v10.5.txt' or
            '10090.protein.actions.v10.5.txt' into a networkx digraph that can be used as the background network for 
            the rest of our functions. IT ONLY KEEPS ACTIVATING AND INHIBITING EDGES FROM THE STRING NETWORK.
            **Intended for Upstream Regulator Analysis functions.

            Args:
                filename: String, the filepath to the STRING database file.
                confidence_filter: A number between 0 and 1000, all interactions with confidece less than this number will be filtered out
        """

        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        DG_universe = create_graph.load_STRING_to_digraph(filename, None, confidence_filter, self.gene_type, self.species)
        self.DG_universe = DG_universe
        
    def load_ndex_from_server(self, UUID, relabel_node_field = None):
    
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
        """
    
        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return
    
        DG_universe = create_graph.load_ndex_from_server(UUID, relabel_node_field, None)
        self.DG_universe = DG_universe
        
    def load_STRING_links(self, filename, confidence_filter = 700):
    
        """
            This function loads a subset of the STRING database from either '9606.protein.links.v10.5.txt' or
            '10090.protein.links.v10.5.txt' into a networkx graph that can be used as the background network for 
            the rest of our functions. **Intended for Heat Propogation and Clustering Analysis functions.

            Args:
                filename: String, the filepath to the STRING database file
                confidence_filter: A number between 0 and 1000, all interactions with confidece less than this number will be filtered out
                species: String, specify either 'human' or 'mouse', should match the specs of your differetial expression database
                translate_to: String, specify either 'symbol' or 'entrez', should match the specs of your differetial expression database
        """
    
        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return
    
        # load STRING
        self.DG_universe = create_graph.load_STRING_links(filename, confidence_filter, self.species, self.gene_type)

        
# --------------------------------------------- LOCALIZATION ------------------------------------------------------#
    
    def localization(self, num_reps = 10, sample_frac = 0.8, method = 'numedges', plot = True, print_counter = False):
    
        """
            Function to calculate localization of your DEGs on your background network.
            Option to compute number of edges (method = 'numedges') or largest connected component (method = 'LLC') 
            localization analysis. Calculates by sampling sub-sections of the focal genes/random set. Percentage to sample
            is set by sample_frac. Option to plot the distributions of random and focal gene localizaiton.
            
            Args:
                num_reps: Int, number of times to randomly sample
                sample_frac: Float, percent of sampled genes
                method: String, to decide which type of localization analysis to run. Options: 'numedges', 'LLC', or 'both'.
                plot: Bool, whether to plot the distributions in the output jupyter notebook cell
                print_counter: Bool, whether to print a counter that tells you which iteration you are on (every 25 iterations).
                               Useful when the num_reps is very high.
                
            Returns: 
                numedges_list: List, the number of edges calculated for each rep, sampling over focal genes. 
                    Empty if method = 'LLC'. 
                numedges_rand: List, the number of edges calculated for each rep, sampling over random genes of 
                    similar degree in the background network. Empty if method = 'LLC'.
                LCC_list: List, the size of the largest connected component, calculated for each rep, sampling over focal genes. 
                    Empty if method = 'numedges'. 
                LCC_rand: List, the size of the largest connected component, calculated for each rep, sampling over random genes of 
                    similar degree in the background network. Empty if method = 'numedges'. 
        """
    
        # make sure user has run all prerequisites
        for item in ['DG_universe', 'DEG_list']:
            if self.check_exists(item) == False:
                return
    
        return stat_analysis.localization(self.DG_universe, self.DEG_list, num_reps, sample_frac, method, plot, print_counter)
        
    def localization_full(self, num_reps = 200, 
                          method = 'LCC', 
                          print_counter = False, 
                          label = 'focal genes',
                          line_height = 0.1,
                          legend_loc = 'upper left'):
    
        """
            Function to calculate localization of an input set of genes (focal_genes) on a background network (Gint).
            Option to compute number of edges (method = 'numedges') or largest connected component (method = 'LLC') 
            localization analysis. DOes no sub-sampling. Plots the distribution of random gene localizaiton, and 
            marks the focal set localization on distribution. Includes p-value of focal set localization.
            
            Args:
                num_reps: Int, number of times to randomly sample
                method: String, to decide which type of localization analysis to run. Options: 'numedges', 'LLC', or 'both'.
                print_counter: Bool, whether to print a counter that tells you which iteration you are on (every 25 iterations).
                               Useful when the num_reps is very high.
                label: String, label for focal genes in graph legend
                line_height: Float, the height of the red line that marks the focal gene localization
                legend_loc: String, relative position of legend in graph. Something similar to 'upper left'.
                
            Returns: 
                numedges_list: List, the number of edges calculated for each rep, over focal genes. 
                    Empty if method = 'LLC'. 
                numedges_rand: List, the number of edges calculated for each rep, over random genes of 
                    similar degree in the background network. Empty if method = 'LLC'.
                LCC_list: List, the size of the largest connected component, calculated for each rep, over focal genes. 
                    Empty if method = 'numedges'. 
                LCC_rand: List, the size of the largest connected component, calculated for each rep, over random genes of 
                    similar degree in the background network. Empty if method = 'numedges'. 
        """
    
        # make sure user has run all prerequisites
        for item in ['DG_universe', 'DEG_list']:
            if self.check_exists(item) == False:
                return
    
        return stat_analysis.localization_full(self.DG_universe, self.DEG_list, num_reps, method, print_counter, label, line_height, legend_loc)

        
        
#------------------------- Heat Propagation --------------------------------#

    def normalized_adj_matrix(self):
    
        """
            This function returns a normalized adjacency matrix.

            Inputs:
                G: NetworkX graph from which to calculate normalized adjacency matrix
                conserve_heat:
                    True: Heat will be conserved (sum of heat vector = 1).  Graph asymmetric
                    False:  Heat will not be conserved.  Graph symmetric.
        """
    
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
                        
        """
            Implements and displays the network propagation for a given graph and seed
            nodes. Additional kwargs are passed to visualizations and visJS_module.

            Inputs:
                num_nodes: Int, the number of the hottest nodes to graph
                edges_width: Int, width of edges in visualized network
                node_size_multiplier: Int, number used to scale the size of all nodes inthe graph
                largest_connected_component: Boolean, whether or not to display largest_connected_component.
                physics_enabled: Boolean, True enables the physics simulation
                node_font_size: Int, font size used to display the gene names of the seed nodes.
                
                ** See visJS2jupyter.visualizations.draw_heat_prop() and visJS2jupyter.visJS_module.visjs_network() for further
                parameter options.

            Returns: VisJS html network plot (iframe) of the heat propagation.
        """
    
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
               
        """
            Creates a visJS2jupyter interactive network in a Jupyter notebook cell that separated the input graph into
            clusters, and colors those clusters either to indicated which cluster they are in, or by log fold change.
            
            Args:
                rad_positions: Boolean, True to separate nodes by cluster, False to leave intermixed
                k: Float, parameter given to networkx.spring_layout() function. If in doubt, leave as None.
                largest_connected_component: Boolean, True to visualize only the largest connected component, False to visualize all nodes.
                alpha: Not currently in use. Functionality may be added later.
                num_its: Not currently in use. Functionality may be added later. 
                num_top_genes: Int, number of genes to visualize 
                cluster_size_cut_off: Int, colors clusters below this size grey.
                remove_stray_nodes: Int, remove clusters of size below cluster_size_cut_off from the visualization.
                r: Float, radius of cluster separation (rad_positions must be True) 
                x_offset: Int, helper that moves clusters around if some are laying on top of each other
                y_offset: Int, helper that moves clusters around if some are laying on top of each other
                node_spacing: Int, increase if there is a lot of overlap between nodes (will happen if there are many nodes in the graph)
                node_size_multiplier: Int, scales the size of each node by this number
                physics_enabled: Boolean, allow interactive, movable nodes and edges. rad_positions must be False.
                node_font_size: Int, font size of the labels marking each seed node
                graph_id: Int, change between visJS calls if you want to visualize multiple networkx int the same jupyter notebook
      
            Returns: The network that will be visualized in the jupyter notebook cell.
        """

        # make sure user has run all prerequisites
        for item in ['DEG_list', 'DG_universe', 'Wprime']:
            if self.check_exists(item) == False:
                return
                
        seed_nodes = [n for n in self.DEG_list if n in self.DG_universe]

        self.node_to_cluster, to_return = heat_and_cluster.draw_clustering(self.DG_universe, seed_nodes,
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
                    node_to_lfc = self.node_to_lfc,
                    **kwargs
                    )
        return to_return
        
        
    def cluster_legend(self, cluster_size_cut_off = 5):
    
        """
            This draws a matplotlib legend displaying the colors associated with each cluster

            Args:
                cluster_size_cut_off: Int, colors clusters below this size grey.

            Returns: A matplotlib figure
        """
    
        # make sure user has run all prerequisites
        for item in ['node_to_cluster']:
            if self.check_exists(item) == False:
                return
        
        node_to_color = heat_and_cluster.assign_colors_to_clusters(self.node_to_cluster, cluster_size_cut_off) # get colors
        node_to_color = pd.Series(node_to_color).drop_duplicates() # clean up color dict
        node_to_cluster = pd.Series(self.node_to_cluster) # clean up cluster dict

        # create a mapping from cluster id to color
        cluster_to_color = {}
        for gene in node_to_color.index.tolist():
            color = node_to_color.loc[gene]
            cluster_id = node_to_cluster.loc[gene]
            cluster_to_color[cluster_id] = color

        # plot the legend
        plt.figure(figsize = (12,6))
        dcount = -1
        for cluster in node_to_cluster.value_counts().keys():
            dcount += 1
            plt.plot([0.1], [dcount], 'o', color = cluster_to_color[cluster])
            plt.annotate('cluster ' + str(cluster), [0.2, dcount - 0.1])

        plt.xlim([0, 3])
        plt.ylim([-5, len(cluster_to_color) + 5])
        
    def draw_legend(self, vmin, vmax, cmap = mpl.cm.bwr, label = 'Units'):
    
        """
            This draws a colormap in a cell that is colored with cmap and ranges from vmin to vmax. 

            Args:
                vmin: Int, the minimum value to display in the colormap
                vmax: Int, the maximum value to display in the colormap
                cmap: matplotlib colormap, used to color the colorbar
                label: String, the text to display below the colorbar

            Returns: A matplotlib colorbar
        """
    
        heat_and_cluster.draw_legend(vmin, vmax, cmap, label)
                    
#------------------ SAVE DATA TO FILE -----------------------------------#

    # write cluster id, seed node Y/N, differential expression log-fold-change, and differential 
    # expression adjusted p-value to a file
    def write_cluster_table(self, to_write_filename):
    
        """
            Writes all of the information we ahve collected in this analysis (gene, cluster id, seed nodes Y/N, log fold change, and p-value).

            Args:
                to_write_filename: String, the file to write the table to.
        """
    
        # make sure user has run all prerequisites
        for item in ['DG_universe', 'DEG_list', 'node_to_cluster', 'node_to_lfc', 'node_to_pvalue']:
            if self.check_exists(item) == False:
                return

        f = open(to_write_filename,"w+")
        f.write('gene,cluster,seed-node,lfc,p-value\n')

        # map nodes to bool value of whether they are a DEG or not
        node_to_DEG_bool = {node:True if node in self.DEG_list else False for node in self.DG_universe}

        for node in self.DG_universe:
            try:
                cluster = self.node_to_cluster[node]
                seed_node = node_to_DEG_bool[node]
                lfc = self.node_to_lfc[node]
                p = self.node_to_pvalue[node]
                
                to_write_string = str(node) + ',' + str(cluster) + ',' + str(seed_node) + ',' + str(lfc) + ',' + str(p) + '\n'
                f.write(to_write_string) 
            except:   
                pass


#------------------ HELPER FUNCTIONS ------------------------------------#  

 
    # item must be string version
    def check_exists(self, item):
    
        """
            Helper functions used to help guide the user as to which order they should call the analysis functions.
            Checks whether the named instance variable has been set yet.

            Args:
                item: String, name of the instance variable whose existance we are checking.
        """

        # re-map it so it stays up to date
        self.string_to_item = {}
        self.string_to_item['gene_type'] = self.gene_type
        self.string_to_item['species'] = self.species
        self.string_to_item['DEG_list'] = self.DEG_list
        self.string_to_item['node_to_lfc'] = self.node_to_lfc
        self.string_to_item['node_to_pvalue'] = self.node_to_lfc
        self.string_to_item['DG_universe'] = self.DG_universe
        self.string_to_item['Wprime'] = self.Wprime
        self.string_to_item['node_to_cluster'] = self.node_to_cluster

        try:
            if (type(self.string_to_item[item]) == type(None)):
                print(self.item_to_message[item])
                return False
        except:
            print('The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
            + '- gene_type\n' \
            + '- species\n' \
            + '- DEG_list\n' \
            + '- node_to_lfc\n' \
            + '- DG_universe\n' \
            + '- Wprime\n' \
            + '- node_to_cluster\n\n')
            return False
        return True