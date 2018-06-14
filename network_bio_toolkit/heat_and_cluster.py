"""
-------------------------------------------
Authors: Mikayla Webster (13webstermj@gmail.com)
         Brin Rosenthal (sbrosenthal@ucsd.edu)
Date: 4/6/18
-------------------------------------------
"""

# our modules
import visJS2jupyter.visJS_module as visJS_module # pip install visJS2jupyter
import visJS2jupyter.visualizations as visualizations # pip install visJS2jupyter

# common packages, most likely already installed
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# uncommon packages required for this analysis
import community # pip install python-louvain

#for local testing
#import sys
#sys.path.append('../../../visJS2jupyter/visJS2jupyter')
#import visualizations
#import visJS_module

# ------------------------- NETWORK VISUALIZATION FUNCTIONS -------------------#

def draw_clustering(DG_universe, seed_nodes,
                    color_lfc = False,
                    rad_positions = True,
                    Wprime = None,
                    k = 0.5,
                    largest_connected_component = True,
                    alpha = 0.5,
                    num_its = 20,
                    num_top_genes = 200,
                    cluster_size_cut_off = 5,
                    remove_stray_nodes = True,
                    r = 1.0,
                    x_offset = 2,
                    y_offset = 2,
                    node_spacing = 700,
                    node_size_multiplier = 4,
                    physics_enabled = False,
                    node_font_size = 45,
                    graph_id = 3,
                    color_map = plt.cm.bwr,
                    alpha_lfc = 1.0, 
				    color_vals_transform = None,
				    ceil_val = 10,
                    color_max_frac = 1.0,
				    color_min_frac = 0.0,
                    node_to_lfc = None,
				    vmin = None,
				    vmax = None,
                    **kwargs
                   ):
                   
    """
        Creates a visJS2jupyter interactive network in a Jupyter notebook cell that separated the input graph into
        clusters, and colors those clusters either to indicated which cluster they are in, or by log fold change.
        
        Args:
            DG_universe: Networkx graph, background network for clustering analysis
            
            seed_nodes: List, seed nodes of heat propagation
            
            color_lfc: Boolean, True to color nodes with log fold change, False to color by cluster
            
            rad_positions: Boolean, True to separate nodes by cluster, False to leave intermixed
            
            Wprime: Adjacency matrix for heat propagation. If Wprime = None, Wprime will be calculated 
                within this function
                
            k: FLoat, parameter given to networkx.spring_layout() function. If in doubt, leave as None.
            
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
            
            color_map: matplotlib.cm.*, a networkx colormap used when coloring with log fold change values (color_lfc must be True)
            
            alpha_lfc: input to visJS_module.return_node_to_color(). (color_lfc must be True)
            
            color_vals_transform: input to visJS_module.return_node_to_color(). (color_lfc must be True)
            
            ceil_val: input to visJS_module.return_node_to_color(). (color_lfc must be True)
            
            color_max_frac: input to visJS_module.return_node_to_color(). (color_lfc must be True)
            
            color_min_frac: input to visJS_module.return_node_to_color(). (color_lfc must be True)
            
            node_to_lfc: Dict, mapping gene names to log fold change values. input to visJS_module.return_node_to_color(). (color_lfc must be True)
            
            vmin: input to visJS_module.return_node_to_color(). (color_lfc must be True)
            
            vmax: input to visJS_module.return_node_to_color(). (color_lfc must be True)

        Returns: The network that will be visualized in the jupyter notebook cell.

    """
    
    # find hottest genes
    if Wprime is None:
        Wprime = visualizations.normalized_adj_matrix(DG_universe)
    
 #   Fnew = visualizations.network_propagation(DG_universe, Wprime, seed_nodes, alpha = alpha, num_its = num_its)
    Fnew = visualizations.network_propagation(nx.Graph(DG_universe), Wprime, seed_nodes)
    top_genes = Fnew.sort_values(ascending = False)[0:num_top_genes].index
    G_top_genes = nx.Graph(DG_universe).subgraph(top_genes) # casting to Graph to match heat prop
    
    # keep only the largest connected component
    G_top_genes = nx.Graph(G_top_genes)
    if largest_connected_component:
        G_top_genes = max(nx.connected_component_subgraphs(G_top_genes), key=len)

    # cluster hottest genes
    node_to_cluster = community.best_partition(G_top_genes)

    nodes = G_top_genes.nodes()
    edges = G_top_genes.edges()
    
    # qucik args check
    if node_to_lfc == None:
        color_lfc = False

    # color based on fold change
    if color_lfc == True:
        # define node colors
        node_to_fld = {n: node_to_lfc[n] for n in nodes if n in node_to_lfc} # keep only those in graph G
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

    # position based on cluster
    if k is None:
        pos = nx.spring_layout(G_top_genes)
    else:
        pos = nx.spring_layout(G_top_genes,k=k)
    
    if rad_positions == True:
        bias_position_by_partition(pos, node_to_cluster, r = r, x_offset = x_offset, y_offset = y_offset); # modifies pos in place

    # set the title of each node
    nx.set_node_attributes(G_top_genes, 'cluster', node_to_cluster)
    node_titles = [str(node[0]) + '<br/>cluster = ' + str(node[1]['cluster'])
           for node in G_top_genes.nodes(data=True)]
    node_titles = dict(zip(nodes, node_titles))

    # label only seed nodes
    node_labels = {n:str(n) if n in seed_nodes else '' for n in G_top_genes.nodes()}

    # change shape of seed nodes
    node_shape = {n:'triangle' if n in seed_nodes else 'dot' for n in G_top_genes.nodes()}

    # create the nodes_dict with all relevant fields
    nodes_dict = [{'id':str(n),
           'color':node_to_color[n],
           'node_shape':node_shape[n],
           'node_label':node_labels[n],
           'title':node_titles[n],
           'x':pos[n][0]*node_spacing,
           'y':pos[n][1]*node_spacing} for n in nodes]

    # map nodes to indices for source/target in edges
    node_map = dict(zip(nodes,range(len(nodes))))

    # create the edges_dict with all relevant fields
    edges_dict = [{'source':node_map[edges[i][0]],
           'target':node_map[edges[i][1]],
           'color':'grey'} for i in range(len(edges))]

    return node_to_cluster, visJS_module.visjs_network(nodes_dict,edges_dict,
                                      node_label_field = 'node_label',
                                      node_size_multiplier = node_size_multiplier,
                                      physics_enabled = physics_enabled,
                                      node_font_size = node_font_size,
                                      graph_id = graph_id,
                                      **kwargs) 

                                      
def draw_legend(vmin, vmax, cmap = mpl.cm.bwr, label = 'Units'):

    """
        This draws a colormap in a cell that is colored with cmap and ranges from vmin to vmax. 

        Args:
            vmin: Int, the minimum value to display in the colormap
            vmax: Int, the maximum value to display in the colormap
            cmap: matplotlib colormap, used to color the colorbar
            label: String, the text to display below the colorbar

        Returns: A matplotlib colorbar

    """

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_axes([0.05, 0.80, 0.9, 0.15])

    # Set the colormap and norm to correspond to vmin and vmax of your original graph
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                    norm=norm,
                                    orientation='horizontal')
    # write label/units below color map
    cb.set_label(label)                                      
                                      
                                      
# ------------------- HELPER FUNCTIONS -----------------------------------------#
                                      
def bias_position_by_partition(pos, partition, r=1.0, x_offset = 2, y_offset = 2):

    '''
        Bias the positions by partition membership, to group them together
    '''
    
    partition = pd.Series(partition)
    
    for p in list(np.unique(partition)):

        focal_nodes = partition[partition == p].index.tolist()
        for n in focal_nodes:
            pos[n][0] += pol2cart(r,float(p)/(np.max(partition) + 1) * x_offset * np.pi)[0]
            pos[n][1] += pol2cart(r,float(p)/(np.max(partition) + 1) * y_offset * np.pi)[1]
            
    return pos

def pol2cart(rho, phi):

    '''
        Helper function for bias_position_by_partition().
    '''
    
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def cart2pol(x, y):

    '''
        Was once a helper function for bias_position_by_partition().
    '''
    
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def assign_colors_to_clusters(node_to_cluster, cluster_size_cut_off = 20, color_list = None):

    """
        Assigns each cluster a color.

        Args:
            node_to_cluster: Dict, output of louvain community() method. Maps genes to their cluster id.
            cluster_size_cut_off: Int, Color clusters with less than this number of genes grey.
            color_list: List, if you would like to supply your own colors rather than the ones we automatically generate
                specify of a list of hexidecimal colors here. When in doubt, leave as None.

        Returns: A dict that maps genes to their color.

    """

    node_to_color = {}
    
    if color_list == None:
        color_list = ['#004033', '#ff0044', '#a300cc', '#3000b3', '#000c59', '#00d9ca', 
         '#00f281', '#00660e', '#88ff00', '#a66f00', '#593000', '#cc5200', '#f22000', '#ff408c', '#731d3f', '#571a66', '#1d3f73', '#40bfff', '#21330d', '#e5c339', '#d96cd2', '#9173e6', '#aacc66', '#736f39', '#403120', '#e5a173', '#e57373', '#d9a3ce', '#40303f', '#acb4e6', '#7c92a6', '#bff2ff', '#b6f2d6', '#d9c7a3', '#806060',
        '#a60016', '#330022', '#004759', '#00735c', '#59000c']
   
    cluster_to_nodes = invert_dict(node_to_cluster)
    
    cluster_to_color = {}
    color_index = 0
    for cluster_id, node_list in cluster_to_nodes.items():
        if len(node_list) < cluster_size_cut_off:
            cluster_to_color[cluster_id] = 'grey'
        else:
            cluster_to_color[cluster_id] = color_list[color_index]
            if color_index < len(color_list)-1:
                color_index = color_index + 1

    for node, cluster_id in node_to_cluster.items():
        node_to_color[node] = cluster_to_color[cluster_id]
    
    return node_to_color

def invert_dict(old_dict):

    """
        Helper function for assign_colors_to_clusters().
    """
    
    inv_dict = {}
    for k, v in old_dict.items():
        inv_dict [v] = inv_dict.get(v, [])
        inv_dict [v].append(k)
    return inv_dict 