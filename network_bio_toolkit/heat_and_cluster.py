"""
-------------------------------------------
Author: Mikayla Webster (13webstermj@gmail.com)
Date: 6/4/18
-------------------------------------------
"""

import community # pip install python-louvain
import networkx as nx
import visJS2jupyter.visJS_module as visJS_module # pip install visJS2jupyter
import visJS2jupyter.visualizations as visualizations # pip install visJS2jupyter
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def draw_clustering(G_DEG, DG_universe, seed_nodes,
                    rad_positions = True,
                    Wprime = None,
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
    
    # find hottest genes
    if Wprime is None:
        Wprime = visualizations.normalized_adj_matrix(G_DEG)
    
    Fnew = visualizations.network_propagation(G_DEG, Wprime, seed_nodes, alpha = alpha, num_its = num_its)
    top_genes = Fnew.sort_values(ascending=False)[0:num_top_genes].index
    G_top_genes = DG_universe.subgraph(top_genes)
    
    # keep only the largest connected component
    G_top_genes = nx.Graph(G_top_genes)
    if largest_connected_component:
        G_top_genes = max(nx.connected_component_subgraphs(G_top_genes), key=len)

    # cluster hottest genes
    node_to_cluster = community.best_partition(G_top_genes)

    nodes = G_top_genes.nodes()
    edges = G_top_genes.edges()

    # color based on cluster
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

    return visJS_module.visjs_network(nodes_dict,edges_dict,
                                      node_label_field = 'node_label',
                                      node_size_multiplier = node_size_multiplier,
                                      physics_enabled = physics_enabled,
                                      node_font_size = node_font_size,
                                      graph_id = graph_id,
                                      **kwargs) 

def bias_position_by_partition(pos, partition, r=1.0, x_offset = 2, y_offset = 2):
    '''
    Bias the positions by partition membership, to group them together
    
    '''
    
    partition = pd.Series(partition)
    
    for p in list(np.unique(partition)):

        focal_nodes = partition[partition==p].index.tolist()
        for n in focal_nodes:
            pos[n][0]+=pol2cart(r,float(p)/(np.max(partition)+1)*x_offset*np.pi)[0]
            pos[n][1]+=pol2cart(r,float(p)/(np.max(partition)+1)*y_offset*np.pi)[1]
            
    return pos

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def assign_colors_to_clusters(node_to_cluster, cluster_size_cut_off = 20, color_list = None):

    node_to_color = {}
    
    if color_list == None:
        color_list = ['#004033', '#ff0044', '#a300cc', '#3000b3', '#000c59', '#00d9ca', 
         '#00f281', '#00660e', '#88ff00', '#a66f00', '#593000', '#cc5200', '#f22000', '#ff408c', '#731d3f', '#571a66', '#1d3f73', '#40bfff', '#21330d', '#e5c339', '#d96cd2', '#9173e6', '#aacc66', '#736f39', '#403120', '#e5a173', '#e57373', '#d9a3ce', '#40303f', '#acb4e6', '#7c92a6', '#bff2ff', '#b6f2d6', '#d9c7a3', '#806060',
        '#a60016', '#330022', '#004759', '#00735c', '#59000c']
   
    cluster_to_nodes = invert_dict(node_to_cluster)
    
    cluster_to_color = {}
    color_index = 0
    for cluster_id, node_list in cluster_to_nodes.iteritems():
        if len(node_list) < cluster_size_cut_off:
            cluster_to_color[cluster_id] = 'grey'
        else:
            cluster_to_color[cluster_id] = color_list[color_index]
            if color_index < len(color_list)-1:
                color_index = color_index + 1

    for node, cluster_id in node_to_cluster.iteritems():
        node_to_color[node] = cluster_to_color[cluster_id]
    
    return node_to_color

def invert_dict(old_dict):
    inv_dict = {}
    for k, v in old_dict.iteritems():
        inv_dict [v] = inv_dict.get(v, [])
        inv_dict [v].append(k)
    return inv_dict 