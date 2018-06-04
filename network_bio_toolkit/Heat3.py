"""
-------------------------------------------
Author: Mikayla Webster (13webstermj@gmail.com)
Date: 4/6/18
-------------------------------------------
"""

import create_graph
import heat_and_cluster
import visJS2jupyter.visualizations as visualizations # pip install visJS2jupyter
import visJS2jupyter.visJS_module as visJS_module
import community # pip install python-louvain
from gprofiler import GProfiler # pip install gprofiler-official

import importlib
importlib.reload(create_graph)
importlib.reload(heat_and_cluster)

#for local testing
#import sys
#sys.path.append('../../../visJS2jupyter/visJS2jupyter')
#import visualizations
#import visJS_module

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

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
        self.node_to_lfc = None
        self.node_to_pvalue = None
        self.DG_universe = None
        self.Wprime = None
        self.node_to_cluster = None
        self.cluster_to_annotation = None
        
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
        self.string_to_item['cluster_to_annotation'] = self.cluster_to_annotation
        
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
                     + 'Or assign your own using Heat_instance.node_to_pvalue\n'
        self.item_to_message['DG_universe'] = 'No background network currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.load_STRING_to_digraph\n' \
                     + 'Or assign your own using Heat_instance.DG_universe\n'
        self.item_to_message['Wprime'] = 'No adjacency matrix currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.normalized_adj_matrix()\n' \
                     + 'Or assign your own using Heat_instance.Wprime\n'
        self.item_to_message['node_to_cluster'] = 'No cluster information currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.draw_clustering_with_annotation()\n' \
                     + 'Or assign your own using Heat_instance.node_to_cluster\n'
        self.item_to_message['cluster_to_annotation'] = 'No annotation information currently on file. Please run the following method:\n' \
                     + ' - Heat_instance.draw_clustering_with_annotation()\n' \
                     + 'Or assign your own using Heat_instance.cluster_to_annotation\n'

            
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
        elif item == 'node_to_lfc': self.node_to_lfc = value
        elif item == 'node_to_pvalue': self.node_to_lfc = value
        elif item == 'DG_universe': self.DG_universe = value
        elif item == 'Wprime': self.Wprime = value
        elif item == 'node_to_cluster': self.node_to_cluster = value
        elif item == 'cluster_to_annotation': self.cluster_to_annotation = value
        
        else:
            print ('The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
            + '- gene_type\n' \
            + '- species\n' \
            + '- DEG_list\n' \
            + '- node_to_lfc\n' \
            + '- node_to_pvalue\n' \
            + '- DG_universe\n' \
            + '- Wprime\n' \
            + '- node_to_cluster\n' \
            + '- cluster_to_annotation\n\n')


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
        DEG_list, node_to_pvalue, node_to_lfc = create_graph.create_DEG_list(filename, None, p_value_filter, p_value_or_adj,
                fold_change_filter, self.gene_type, gene_column_header, p_value_column_header, fold_change_column_header, sep,
                return_full_values = True)
        self.DEG_list = DEG_list
        self.node_to_lfc = node_to_lfc
        self.node_to_pvalue = node_to_pvalue
        
        
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
                    
#----------------------- ANNOTATION -------------------------------------#
     
          
    def draw_clustering_with_annotation(self, 
                        num_nodes = 500,
                        annotation = True,
                        node_spacing = 700,
                        node_size_multiplier = 5,
                        physics_enabled = False,
                        node_font_size = 20,
                        graph_id = 4,
                        node_size = 7,
                        cluster_size_cutoff = 5,
                        color_lfc = True,
                        remove_stray_nodes = True,
                        color_map = plt.cm.bwr,
                        alpha_lfc = 1.0, 
                        color_vals_transform = None,
                        ceil_val = 10,
                        color_max_frac = 1.0,
                        color_min_frac = 0.0,
                        vmin = None,
                        vmax = None,
                        **kwargs):
                        
        # make sure user has run all prerequisites
        for item in ['DEG_list', 'DG_universe', 'Wprime', 'node_to_lfc']:
            if self.check_exists(item) == False:
                return
        
        # run heat prop
        seed_nodes = [n for n in self.DEG_list if n in self.DG_universe]
        Fnew = visualizations.network_propagation(nx.Graph(self.DG_universe), self.Wprime, seed_nodes)
        top_genes = Fnew.sort_values(ascending=False)[0:num_nodes].index
        G_top_genes = nx.Graph(self.DG_universe).subgraph(top_genes) # casting to Graph to match heat prop

        # keep only the largest connected component
        G_top_genes = nx.Graph(G_top_genes)
        G_top_genes = max(nx.connected_component_subgraphs(G_top_genes), key=len)

        # cluster hottest genes
        node_to_cluster = community.best_partition(G_top_genes)
        
        # annotation
        cluster_to_annotation = self.get_annotations(node_to_cluster, top_genes)
        cluster_to_node = heat_and_cluster.invert_dict(node_to_cluster)
        
        # create edges between each node in a cluster and that cluster's annotation node
        edge_list = []
        for cluster_id, node_group in cluster_to_node.items():
            anno = cluster_to_annotation[cluster_id]
            node_to_cluster[anno] = cluster_id
            anno_list = [anno]*len(node_group)
            new_edges = list(zip(anno_list, node_group))
            edge_list.extend(new_edges) 
            
        # incorporate annotation into graph
        G_top_genes.add_edges_from(edge_list)
        nodes = G_top_genes.nodes()
        edges = G_top_genes.edges()

        # position based on cluster
        pos = nx.spring_layout(G_top_genes)
        heat_and_cluster.bias_position_by_partition(pos, node_to_cluster, r = 0.5); # modifies pos in place

        # color based on fold change
        if color_lfc == True:
            # define node colors
            node_to_fld = {n: self.node_to_lfc[n] for n in nodes if n in self.node_to_lfc} # keep only those in graph G
            nx.set_node_attributes(G_top_genes, 'fold_change', 0) # give all nodes a default fold change of zero
            nx.set_node_attributes(G_top_genes, 'fold_change', node_to_fld) # overwrite with actual fold change for the nodes that have one
            node_to_color = visJS_module.return_node_to_color(G_top_genes, field_to_map = 'fold_change', 
                                                                cmap = color_map, 
                                                                alpha = alpha_lfc, 
                                                                color_vals_transform = color_vals_transform,
                                                                ceil_val = ceil_val,
                                                                color_max_frac = color_max_frac,
                                                                color_min_frac = color_min_frac,
                                                                vmin = vmin,
                                                                vmax = vmax)
        else: # color based on cluster
            node_to_color = assign_colors_to_clusters(node_to_cluster, cluster_size_cut_off)

        # remove stray nodes
        if remove_stray_nodes == True:
            not_grey_list = [n for n in nodes if node_to_color[n] != 'grey']
            node_to_cluster = {n:node_to_cluster[n] for n in nodes if n in not_grey_list}

            G_top_genes = nx.subgraph(G_top_genes, not_grey_list)
            nodes = G_top_genes.nodes()
            edges = G_top_genes.edges()

        # set the annotation nodes to white
        for k,v in cluster_to_annotation.items():
            node_to_color[v] = 'white'
            
        # set the title of each node
        nx.set_node_attributes(G_top_genes, 'cluster', node_to_cluster)
        node_titles = [str(node[0]) + '<br/>cluster = ' + str(node[1]['cluster']) for node in G_top_genes.nodes(data=True)]
        node_titles = dict(zip(nodes, node_titles))
        
        # label only seed nodes
        if annotation == False:
            node_labels = {n:str(n) if n in seed_nodes else '' for n in G_top_genes.nodes()}
        else:
            node_labels = {n:'' for n in G_top_genes.nodes()}
            for k,v in cluster_to_annotation.items():
                node_labels[v] = v
                
        # change shape of seed nodes
        node_shape = {n:'triangle' if n in seed_nodes else 'dot' for n in G_top_genes.nodes()}
        
        # set node size
        node_to_size = {n:node_size for n in G_top_genes.nodes()}
        #for k,v in cluster_to_annotation.items():
        #    node_to_size[v] = 1

        # create the nodes_dict with all relevant fields
        nodes_dict = [{'id':str(n),
               'color':node_to_color[n],
               'node_shape':node_shape[n],
               'node_label':node_labels[n],
               'node_size': node_to_size[n],
               'title':node_titles[n],
               'x':pos[n][0]*node_spacing,
               'y':pos[n][1]*node_spacing} for n in nodes]

        # map nodes to indices for source/target in edges
        node_map = dict(zip(nodes,range(len(nodes))))

        # decide edge color
        edge_list = [tuple(reversed(edge)) for edge in edge_list] # tuples are in the wrong direction for some reason
        edge_to_color = {edge:'white' if edge in edge_list else 'grey' for edge in edges} # make edges to annotation nodes white

        # move annotation nodes to front of list (puts annotation edges behind all other edges in the graph)
        nodes = list(reversed(nodes))
        for edge in edge_list:
            edges.insert(0, edges.pop(edges.index(edge)))
            
        # create the edges_dict with all relevant fields
        edges_dict = [{'source':node_map[edges[i][0]],
               'target':node_map[edges[i][1]],
               'color':edge_to_color[edges[i]]} for i in range(len(edges))]
            
        # save our node_to_cluster dict and 
        self.node_to_cluster = node_to_cluster
        self.cluster_to_annotation = cluster_to_annotation

        return visJS_module.visjs_network(nodes_dict,edges_dict,
                                          node_label_field = 'node_label',
                                          node_size_field = 'node_size',
                                          physics_enabled = physics_enabled,
                                          node_font_size = 20,
                                          graph_id = graph_id,
                                          **kwargs)
                                          
    def cluster_legend(self, cluster_size_cut_off = 5):
    
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
                                          
#------------------ SAVE DATA TO FILE -----------------------------------#

    # write cluster id, seed node Y/N, annotation, differential expression log-fold-change, and differential 
    # expression adjusted p-value to a file
    def write_cluster_table(self, to_write_filename):
    
        # make sure user has run all prerequisites
        for item in ['DG_universe', 'DEG_list', 'node_to_cluster', 'node_to_lfc', 'node_to_pvalue', 'cluster_to_annotation']:
            if self.check_exists(item) == False:
                return

        f = open(to_write_filename,"w+")
        f.write('gene,cluster,seed-node,annotation,lfc,p-value\n')

        # map nodes to bool value of whether they are a DEG or not
        node_to_DEG_bool = {node:True if node in self.DEG_list else False for node in self.DG_universe}

        for node in self.DG_universe:
            try:
                cluster = self.node_to_cluster[node]
                seed_node = node_to_DEG_bool[node]
                anno = self.cluster_to_annotation[self.node_to_cluster[node]]
                lfc = self.node_to_lfc[node]
                p = self.node_to_pvalue[node]
                
                to_write_string = str(node) + ',' + str(cluster) + ',' + str(seed_node) + ',' + str(anno) + ',' + str(lfc) + ',' + str(p) + '\n'
                f.write(to_write_string) 
            except:
                pass
        

#------------------ HELPER FUNCTIONS ------------------------------------#  

    def get_annotations(self, node_to_cluster, top_nodes, top = None):
        
        cluster_to_node = heat_and_cluster.invert_dict(node_to_cluster)
        counter = 0
        
        if top == None:
            top = len(cluster_to_node) - 1
        
        gp = GProfiler('cluster_annotator', want_header = True)
        
        cluster_to_annotation = {}
        for cluster_id, group in cluster_to_node.items():
            
            print('Annotating cluster ' + str(cluster_id) + ' of ' + str(top) + '...')
            
            if self.species == 'human':
                org = 'hsapiens'
            else:
                org = 'mmusculus'
            df = pd.DataFrame(gp.gprofile(group, organism = org, custom_bg = top_nodes))
            df.columns = df.loc[0]
            df.drop(0)
            annotation = df['t name'][1]
            
            cluster_to_annotation[cluster_id] = annotation
            
            if counter == top:
                print('Done!')
                return cluster_to_annotation
            counter = counter + 1

    # item must be string version
    def check_exists(self, item):

        # re-map it so it stays up to date
        self.string_to_item = {}
        self.string_to_item['gene_type'] = self.gene_type
        self.string_to_item['species'] = self.species
        self.string_to_item['DEG_list'] = self.DEG_list
        self.string_to_item['node_to_lfc'] = self.node_to_lfc
        self.string_to_item['node_to_pvalue'] = self.node_to_pvalue
        self.string_to_item['DG_universe'] = self.DG_universe
        self.string_to_item['Wprime'] = self.Wprime
        self.string_to_item['node_to_cluster'] = self.node_to_cluster
        self.string_to_item['cluster_to_annotation'] = self.cluster_to_annotation

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
            + '- node_to_pvalue\n' \
            + '- DG_universe\n' \
            + '- Wprime\n' \
            + '- node_to_cluster\n' \
            + '- cluster_to_annotation\n\n')
            return False
        return True