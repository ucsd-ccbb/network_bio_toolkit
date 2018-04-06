from __future__ import print_function
import pandas as pd
import networkx as nx
import numpy as np
import visJS2jupyter.visJS_module as visJS_module
import visJS2jupyter.visualizations as visualizations
import matplotlib as plt
import json
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
#import visJS_module # use this for local testing
import visJS2jupyter.visJS_module as visJS_module
import community # pip install python-louvain
import networkx as nx
import visJS2jupyter.visJS_module as visJS_module
import visJS2jupyter.visualizations as visualizations
import matplotlib as mpl
import random
import operator
import matplotlib.pyplot as plt
import numpy as np
import copy


# ------------------ HEAT PROP ---------------------------------#

def normalized_adj_matrix(G,conserve_heat=True,weighted=False):
    return visualizations.normalized_adj_matrix(G,conserve_heat,weighted)
	
def network_propagation(G,Wprime,seed_nodes,alpha=.5, num_its=20):
    return network_propagation(G,Wprime,seed_nodes,alpha, num_its)

def draw_heat_prop(G, DEG_list, 
                   num_nodes = 200,
                   Wprime = None,
                   k=.7,
                   edge_smooth_enabled=True,
                   edge_smooth_type='bezier',
                   largest_connected_component=True,
                   node_size_multiplier=2,
                   highlight_nodes=None, 
                   **kwargs):
    
    DEG_list = [n for n in DEG_list if n in G]
    G = nx.Graph(G)
	
    if highlight_nodes == None:
	    highlight_nodes = DEG_list
		
    return vis_draw_heat_prop(G, DEG_list, 
                                         num_nodes = num_nodes,
                                         Wprime = Wprime,
                                         k = k,
                                         edge_smooth_enabled = edge_smooth_enabled,
                                         edge_smooth_type = edge_smooth_type,
                                         largest_connected_component = largest_connected_component,
                                         node_size_multiplier = node_size_multiplier,
                                         highlight_nodes = highlight_nodes, 
                                         **kwargs)
										 
										 
def vis_draw_heat_prop(G, seed_nodes, random_walk = True,
                   edge_cmap=plt.cm.autumn_r,
                   export_file='heat_prop.json',
                   export_network=False,
                   highlight_nodes=None,
                   k=None,
                   largest_connected_component=False,
                   node_cmap=plt.cm.autumn_r,
                   node_size=10,
                   num_nodes=None,
                   physics_enabled=False,
                   Wprime=None,
                   **kwargs):
    # check for invalid nodes in seed_nodes
    invalid_nodes = [node for node in seed_nodes if node not in G.nodes()]
    for node in invalid_nodes:
        print ('Node {} not in graph'.format(node))
    if invalid_nodes:
        return

    # perform the network propagation
    if random_walk == True: # perform random walk style heat propagation
        if Wprime is None:
            Wprime = visualizations.normalized_adj_matrix(G)    
        prop_graph = visualizations.network_propagation(G, Wprime, seed_nodes).to_dict()
        nx.set_node_attributes(G, name = 'node_heat', values = prop_graph)
    #else: # perform diffusion style heat propagation
    #    if type(seed_nodes) != dict:
    #        print('When parameter random_walk = False, parameter seed_nodes must be a dict')
    #        return -1
    #    heat_kernel = scipy_heatKernel.SciPYKernel(G) # need a graph
    #    diffused_heats = heat_kernel.diffuse(seed_nodes) # need seed_to_heat mapping
    #    nx.set_node_attributes(G, name = 'node_heat', values = dict(diffused_heats))

    # find top num_nodes hottest nodes and connected component if requested
    G = visualizations.set_num_nodes(G,num_nodes)
    if largest_connected_component:
        G = max(nx.connected_component_subgraphs(G), key=len)
    nodes = list(G.nodes())
    edges = list(G.edges())

    # check for empty nodes and edges after getting subgraph of G
    if not nodes:
        print ('There are no nodes in the graph. Try increasing num_nodes.')
        return
    if not edges:
        print ('There are no edges in the graph. Try increasing num_nodes.')
        return

    # set the position of each node
    if k is None:
        pos = nx.spring_layout(G)
    else:
        pos = nx.spring_layout(G,k=k)

    xpos,ypos=zip(*pos.values())
    nx.set_node_attributes(G, name = 'xpos', values = dict(zip(pos.keys(),[x*1000 for x in xpos])))
    nx.set_node_attributes(G, name = 'ypos', values = dict(zip(pos.keys(),[y*1000 for y in ypos])))

    # set the border width of nodes
    if 'node_border_width' not in kwargs.keys():
        kwargs['node_border_width'] = 2

    border_width = {}
    for n in nodes:
        if n in seed_nodes:
            border_width[n] = kwargs['node_border_width']
        elif highlight_nodes is not None and n in highlight_nodes:
            border_width[n] = kwargs['node_border_width']
        else:
            border_width[n] = 0

    nx.set_node_attributes(G, name = 'nodeOutline', values = border_width)

    # set the shape of each node
    nodes_shape=[]
    for node in G.nodes():
        if node in seed_nodes:
            nodes_shape.append('triangle')
        else:
            nodes_shape.append('dot')
    node_to_shape=dict(zip(G.nodes(),nodes_shape))
    nx.set_node_attributes(G, name = 'nodeShape', values = node_to_shape)

    # add a field for node labels
    if highlight_nodes:
        node_labels = {}
        for node in nodes:
            if node in seed_nodes:
                node_labels[node] = str(node)
            elif node in highlight_nodes:
                node_labels[node] = str(node)
            else:
                node_labels[node] = ''
    else:
        node_labels = {n:str(n) for n in nodes}

    nx.set_node_attributes(G, name = 'nodeLabel', values = node_labels)

    # set title for each node
    node_titles = [str(node[0]) + '<br/>heat = ' + str(round(node[1]['node_heat'],5))
                   for node in G.nodes(data=True)]
    node_titles = dict(zip(G.nodes(),node_titles))
    nx.set_node_attributes(G, name = 'nodeTitle', values = node_titles)

    # set color of each node
    node_to_color = visJS_module.return_node_to_color(G,
                                                      field_to_map='node_heat',
                                                      cmap=node_cmap,
                                                      color_vals_transform='log')

    # set heat value of edge based off hottest connecting node's value
    node_attr = nx.get_node_attributes(G,'node_heat')
    edge_weights = {}
    for e in edges:
        if node_attr[e[0]] > node_attr[e[1]]:
            edge_weights[e] = node_attr[e[0]]
        else:
            edge_weights[e] = node_attr[e[1]]

    nx.set_edge_attributes(G, name = 'edge_weight', values = edge_weights)

    # set color of each edge
    edge_to_color = visJS_module.return_edge_to_color(G,
                                                      field_to_map='edge_weight',
                                                      cmap=edge_cmap,
                                                      color_vals_transform='log')

    # create the nodes_dict with all relevant fields
    nodes_dict = [{'id':str(n),
                   'border_width':border_width[n],
                   'degree':G.degree(n),
                   'color':node_to_color[n],
                   'node_label':node_labels[n],
                   'node_size':node_size,
                   'node_shape':node_to_shape[n],
                   'title':node_titles[n],
                   'x':np.float64(pos[n][0]).item()*1000,
                   'y':np.float64(pos[n][1]).item()*1000} for n in nodes]

    # map nodes to indices for source/target in edges
    node_map = dict(zip(nodes,range(len(nodes))))

    # create the edges_dict with all relevant fields
    edges_dict = [{'source':node_map[edges[i][0]],
                   'target':node_map[edges[i][1]],
                   'color':edge_to_color[edges[i]]} for i in range(len(edges))]

    # set node_size_multiplier to increase node size as graph gets smaller
    if 'node_size_multiplier' not in kwargs.keys():
        if len(nodes) > 500:
            kwargs['node_size_multiplier'] = 3
        elif len(nodes) > 200:
            kwargs['node_size_multiplier'] = 5
        else:
            kwargs['node_size_multiplier'] = 7

    kwargs['physics_enabled'] = physics_enabled

    # if node hovering color not set, set default to black
    if 'node_color_hover_background' not in kwargs.keys():
        kwargs['node_color_hover_background'] = 'black'

    # node size determined by size in nodes_dict, not by id
    if 'node_size_field' not in kwargs.keys():
        kwargs['node_size_field'] = 'node_size'

    # node label determined by value in nodes_dict
    if 'node_label_field' not in kwargs.keys():
        kwargs['node_label_field'] = 'node_label'

    # export the network to JSON for Cytoscape
    if export_network:
        node_colors = map_node_to_color(G,'node_heat',True)
        nx.set_node_attributes(G, name = 'nodeColor', values = node_colors)
        edge_colors = map_edge_to_color(G,'edge_weight',True)
        nx.set_edge_attributes(G, name = 'edgeColor', values = edge_colors)
        visJS_module.export_to_cytoscape(G = G,export_file = export_file)

    return visJS_module.visjs_network(nodes_dict,edges_dict,**kwargs)


# ----------------- CLUSTERING ANNOTATION ---------------------#

def draw_clustering(G, DEG_list,
                   cluster_size_cut_off = 20,
                   num_nodes = 200,
                   Wprime = None,
                   largest_connected_component=True,
                   node_size_multiplier=2,
                   highlight_nodes=None, 
                   **kwargs):
    
    DEG_list = [n for n in DEG_list if n in G]
    G = nx.Graph(G)
	
    if highlight_nodes == None:
	    highlight_nodes = DEG_list
		
    return vis_draw_clustering(G, DEG_list, cluster_size_cut_off,
                                         num_nodes = num_nodes,
                                         Wprime = Wprime,
                                         largest_connected_component = largest_connected_component,
                                         node_size_multiplier = node_size_multiplier,
                                         highlight_nodes = highlight_nodes, 
                                         **kwargs)
										 
	

def vis_draw_clustering(G, seed_nodes, 
                   cluster_size_cut_off = 20,
                   edge_cmap=plt.cm.autumn_r,
                   export_file='heat_prop.json',
                   random_walk = True,
                   export_network=False,
                   highlight_nodes=None,
                   k=None,
                   largest_connected_component=False,
                   node_cmap=plt.cm.autumn_r,
                   node_size=10,
                   num_nodes=None,
                   physics_enabled=False,
                   Wprime=None,
                   **kwargs):
    '''
    Implements and displays the network propagation for a given graph and seed
    nodes. Additional kwargs are passed to visJS_module.

    Inputs:
        - G: a networkX graph
        - seed_nodes: nodes on which to initialize the simulation (must be a dict if random_walk = False)
		- random_walk: True to perform a random walk style heat propagation, False to perform a diffusion style one.
        - edge_cmap: matplotlib colormap for edges, default: matplotlib.cm.autumn_r
        - export_file: JSON file to export graph data, default: 'graph_overlap.json'
        - export_network: export network to Cytoscape, default: False
        - highlight_nodes: list of nodes to place borders around, default: None
        - k: float, optimal distance between nodes for nx.spring_layout(), default: None
        - largest_connected_component: boolean, whether or not to display largest_connected_component,
                                       default: False
        - node_cmap: matplotlib colormap for nodes, default: matplotlib.cm.autumn_r
        - node_size: size of nodes, default: 10
        - num_nodes: the number of the hottest nodes to graph, default: None (all nodes will be graphed)
        - physics_enabled: enable physics simulation, default: False
        - Wprime: normalized adjacency matrix (from function normalized_adj_matrix())

    Returns:
        - VisJS html network plot (iframe) of the heat propagation.
    '''

    # check for invalid nodes in seed_nodes
    invalid_nodes = [node for node in seed_nodes if node not in G.nodes()]
    for node in invalid_nodes:
        print ('Node {} not in graph'.format(node))
    if invalid_nodes:
        return

    # perform the network propagation
    if random_walk == True: # perform random walk style heat propagation
        if Wprime is None:
            Wprime = visualizations.normalized_adj_matrix(G)    
        prop_graph = visualizations.network_propagation(G, Wprime, seed_nodes).to_dict()
        nx.set_node_attributes(G, name = 'node_heat', values = prop_graph)
    #else: # perform diffusion style heat propagation
    #    if type(seed_nodes) != dict:
    #        print('When parameter random_walk = False, parameter seed_nodes must be a dict')
    #        return -1
    #    heat_kernel = scipy_heatKernel.SciPYKernel(G) # need a graph
    #    diffused_heats = heat_kernel.diffuse(seed_nodes) # need seed_to_heat mapping
    #    nx.set_node_attributes(G, name = 'node_heat', values = dict(diffused_heats))

    # find top num_nodes hottest nodes and connected component if requested
    G = visualizations.set_num_nodes(G,num_nodes)
    if largest_connected_component:
        G = max(nx.connected_component_subgraphs(G), key=len)
    nodes = list(G.nodes())
    edges = list(G.edges())

    # check for empty nodes and edges after getting subgraph of G
    if not nodes:
        print ('There are no nodes in the graph. Try increasing num_nodes.')
        return
    if not edges:
        print ('There are no edges in the graph. Try increasing num_nodes.')
        return

    # set the position of each node
    if k is None:
        pos = nx.spring_layout(G)
    else:
        pos = nx.spring_layout(G,k=k)

    xpos,ypos=zip(*pos.values())
    nx.set_node_attributes(G, name = 'xpos', values = dict(zip(pos.keys(),[x*1000 for x in xpos])))
    nx.set_node_attributes(G, name = 'ypos', values = dict(zip(pos.keys(),[y*1000 for y in ypos])))

    # set the border width of nodes
    if 'node_border_width' not in kwargs.keys():
        kwargs['node_border_width'] = 2

    border_width = {}
    for n in nodes:
        if n in seed_nodes:
            border_width[n] = kwargs['node_border_width']
        elif highlight_nodes is not None and n in highlight_nodes:
            border_width[n] = kwargs['node_border_width']
        else:
            border_width[n] = 0

    nx.set_node_attributes(G, name = 'nodeOutline', values = border_width)

    # set the shape of each node
    nodes_shape=[]
    for node in G.nodes():
        if node in seed_nodes:
            nodes_shape.append('triangle')
        else:
            nodes_shape.append('dot')
    node_to_shape=dict(zip(G.nodes(),nodes_shape))
    nx.set_node_attributes(G, name = 'nodeShape', values = node_to_shape)

    # add a field for node labels
    if highlight_nodes:
        node_labels = {}
        for node in nodes:
            if node in seed_nodes:
                node_labels[node] = str(node)
            elif node in highlight_nodes:
                node_labels[node] = str(node)
            else:
                node_labels[node] = ''
    else:
        node_labels = {n:str(n) for n in nodes}

    nx.set_node_attributes(G, name = 'nodeLabel', values = node_labels)

    # set color of each node ####################################################################
    # G = nx.Graph(G)
    node_to_cluster = community.best_partition(G)
    nx.set_node_attributes(G, name = 'cluster', values = node_to_cluster)
    node_to_color = {}
    node_to_color = assign_colors_to_clusters(node_to_cluster, cluster_size_cut_off = cluster_size_cut_off, color_list = None)
    
    # set the title of each node
    node_titles = [str(node[0]) + '<br/>cluster = ' + str(node[1]['cluster'])
                   for node in G.nodes(data=True)]
    node_titles = dict(zip(nodes,node_titles))
    nx.set_node_attributes(G, name = 'nodeTitle', values = node_titles)

    # create the nodes_dict with all relevant fields
    nodes_dict = [{'id':str(n),
 #                  'border_width':border_width[n],
                   'degree':G.degree(n),
                   'color':node_to_color[n],
                   'node_label':node_labels[n],
                   'node_size':node_size,
                   'node_shape':node_to_shape[n],
                   'title':node_titles[n],
                   'x':np.float64(pos[n][0]).item()*1000,
                   'y':np.float64(pos[n][1]).item()*1000} for n in nodes]

    # map nodes to indices for source/target in edges
    node_map = dict(zip(nodes,range(len(nodes))))

    # create the edges_dict with all relevant fields
    edges_dict = [{'source':node_map[edges[i][0]],
                   'target':node_map[edges[i][1]],
                   'color':'grey'} for i in range(len(edges))]

    # set node_size_multiplier to increase node size as graph gets smaller
    if 'node_size_multiplier' not in kwargs.keys():
        if len(nodes) > 500:
            kwargs['node_size_multiplier'] = 3
        elif len(nodes) > 200:
            kwargs['node_size_multiplier'] = 5
        else:
            kwargs['node_size_multiplier'] = 7

    kwargs['physics_enabled'] = physics_enabled

    # if node hovering color not set, set default to black
    if 'node_color_hover_background' not in kwargs.keys():
        kwargs['node_color_hover_background'] = 'black'

    # node size determined by size in nodes_dict, not by id
    if 'node_size_field' not in kwargs.keys():
        kwargs['node_size_field'] = 'node_size'

    # node label determined by value in nodes_dict
    if 'node_label_field' not in kwargs.keys():
        kwargs['node_label_field'] = 'node_label'

    # export the network to JSON for Cytoscape
    if export_network:
        node_colors = map_node_to_color(G,'node_heat',True)
        nx.set_node_attributes(G, name = 'nodeColor', values = node_colors)
        edge_colors = map_edge_to_color(G,'edge_weight',True)
        nx.set_edge_attributes(G, name = 'edgeColor', values = edge_colors)
        visJS_module.export_to_cytoscape(G = G,export_file = export_file)

    return visJS_module.visjs_network(nodes_dict,edges_dict,**kwargs)
	
	
def invert_dict(old_dict):
    inv_dict = {}
    for k, v in old_dict.iteritems():
        inv_dict [v] = inv_dict.get(v, [])
        inv_dict [v].append(k)
    return inv_dict 
	
def assign_colors_to_clusters(node_to_cluster, cluster_size_cut_off = 20, color_list = None):

    node_to_color = {}
    
    if color_list == None:
       # color_list = ["#00FF80", "#FF00FF", "#FF8000", "#FF3300", "#FFD11A", "#CC00CC", "#00FFFF"]
	   color_list = ['#ff0044', '#a300cc', '#3000b3', '#000c59', '#00d9ca', 
	   '#004033', '#00f281', '#00660e', '#88ff00', '#a66f00', '#593000', '#cc5200', '#f22000', '#ff408c', '#731d3f', '#571a66', '#1d3f73', '#40bfff', '#21330d', '#e5c339', '#d96cd2', '#9173e6', '#aacc66', '#736f39', '#403120', '#e5a173', '#e57373', '#d9a3ce', '#40303f', '#acb4e6', '#7c92a6', '#bff2ff', '#b6f2d6', '#d9c7a3', '#806060',
	   '#a60016', '#330022', '#004759', '#00735c', '#59000c']
	   
	   # color_list = [#ffbfbf, #ff9f40, #e5e600, #33cc80, #266099, #7b698c, #f23d97, #997373, #8c5823, #4c4d00, #ace6e6, #102840, #502080, #330d20, #d93636, #4c3013, #8fe639, #134d4d, #b6b6f2, #d96cd9, #661a1a, #808060, #408000, #00e6e6, #8080ff, #ff00ff, #f2d4b6, #bfbf60, #8fbf8f, #00cccc, #3939e6, #ffbfdf, #66594d, #ffff00, #104010, #409fff, #13134d, #66334d]
   
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



# ----------------- RADIAL SPACING ANNOTATION ----------------- #

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def bias_position_by_partition(pos, partition, r=0.5):
    '''
    Bias the positions by partition membership, to group them together

    '''
    return_pos = copy.deepcopy(pos)
    partition = pd.Series(partition)

    for p in list(np.unique(partition)):

        focal_nodes = partition[ partition ==p].index.tolist()
        for n in focal_nodes:
            to_add = pol2cart(r ,float(p ) /(np.max(partition ) +1 ) * 2 *np.pi)
            return_pos[n][0 ] = pos[n][0 ] + to_add[0]
            return_pos[n][1 ] = pos[n][1 ] + to_add[1]

    return return_pos
	
	
def draw_rad_clustering(G, DEG_list,
                   cluster_size_cut_off = 20,
                   node_spacing = 3000,
                   r = 0.5,
                   num_nodes = 200,
                   Wprime = None,
                   largest_connected_component=True,
                   node_size_multiplier=2,
                   highlight_nodes=None, 
                   **kwargs):
    
    DEG_list = [n for n in DEG_list if n in G]
    G = nx.Graph(G)
	
    if highlight_nodes == None:
	    highlight_nodes = DEG_list
		
    return vis_draw_rad_clustering(G, DEG_list, cluster_size_cut_off,
                                         node_spacing = node_spacing,
                                         r = r,
                                         num_nodes = num_nodes,
                                         Wprime = Wprime,
                                         largest_connected_component = largest_connected_component,
                                         node_size_multiplier = node_size_multiplier,
                                         highlight_nodes = highlight_nodes, 
                                         **kwargs)
	
	
def vis_draw_rad_clustering(G, seed_nodes,
                   cluster_size_cut_off = 20, 
                   node_spacing = 3000,
                   r = 0.5,
                   num_nodes=None,
                   edge_cmap=plt.cm.autumn_r,
                   export_file='heat_prop.json',
                   random_walk = True,
                   export_network=False,
                   highlight_nodes=None,
                   k=None,
                   largest_connected_component=False,
                   node_cmap=plt.cm.autumn_r,
                   node_size=10,
                   physics_enabled=False,
                   Wprime=None,
                   **kwargs):
    '''
    Implements and displays the network propagation for a given graph and seed
    nodes. Additional kwargs are passed to visJS_module.

    Inputs:
        - G: a networkX graph
        - seed_nodes: nodes on which to initialize the simulation (must be a dict if random_walk = False)
		- random_walk: True to perform a random walk style heat propagation, False to perform a diffusion style one.
        - edge_cmap: matplotlib colormap for edges, default: matplotlib.cm.autumn_r
        - export_file: JSON file to export graph data, default: 'graph_overlap.json'
        - export_network: export network to Cytoscape, default: False
        - highlight_nodes: list of nodes to place borders around, default: None
        - k: float, optimal distance between nodes for nx.spring_layout(), default: None
        - largest_connected_component: boolean, whether or not to display largest_connected_component,
                                       default: False
        - node_cmap: matplotlib colormap for nodes, default: matplotlib.cm.autumn_r
        - node_size: size of nodes, default: 10
        - num_nodes: the number of the hottest nodes to graph, default: None (all nodes will be graphed)
        - physics_enabled: enable physics simulation, default: False
        - Wprime: normalized adjacency matrix (from function normalized_adj_matrix())

    Returns:
        - VisJS html network plot (iframe) of the heat propagation.
    '''

    # check for invalid nodes in seed_nodes
    invalid_nodes = [node for node in seed_nodes if node not in G.nodes()]
    for node in invalid_nodes:
        print ('Node {} not in graph'.format(node))
    if invalid_nodes:
        return

    # perform the network propagation
    if random_walk == True: # perform random walk style heat propagation
        if Wprime is None:
            Wprime = visualizations.normalized_adj_matrix(G)    
        prop_graph = visualizations.network_propagation(G, Wprime, seed_nodes).to_dict()
        nx.set_node_attributes(G, name = 'node_heat', values = prop_graph)
    #else: # perform diffusion style heat propagation
    #    if type(seed_nodes) != dict:
    #        print('When parameter random_walk = False, parameter seed_nodes must be a dict')
    #        return -1
    #    heat_kernel = scipy_heatKernel.SciPYKernel(G) # need a graph
    #    diffused_heats = heat_kernel.diffuse(seed_nodes) # need seed_to_heat mapping
    #    nx.set_node_attributes(G, name = 'node_heat', values = dict(diffused_heats))

    # find top num_nodes hottest nodes and connected component if requested
    G = visualizations.set_num_nodes(G,num_nodes)
    if largest_connected_component:
        G = max(nx.connected_component_subgraphs(G), key=len)
    nodes = list(G.nodes())
    edges = list(G.edges())

    # check for empty nodes and edges after getting subgraph of G
    if not nodes:
        print ('There are no nodes in the graph. Try increasing num_nodes.')
        return
    if not edges:
        print ('There are no edges in the graph. Try increasing num_nodes.')
        return

 #   # set the border width of nodes
 #   if 'node_border_width' not in kwargs.keys():
 #       kwargs['node_border_width'] = 2

 #   border_width = {}
 #   for n in nodes:
 #       if n in seed_nodes:
 #           border_width[n] = kwargs['node_border_width']
 #       elif highlight_nodes is not None and n in highlight_nodes:
 #           border_width[n] = kwargs['node_border_width']
 #       else:
 #           border_width[n] = 0

 #   nx.set_node_attributes(G, name = 'nodeOutline', values = border_width)

    # set the shape of each node
    nodes_shape=[]
    for node in G.nodes():
        if node in seed_nodes:
            nodes_shape.append('triangle')
        else:
            nodes_shape.append('dot')
    node_to_shape=dict(zip(G.nodes(),nodes_shape))
    nx.set_node_attributes(G, name = 'nodeShape', values = node_to_shape)

    # add a field for node labels
    if highlight_nodes:
        node_labels = {}
        for node in nodes:
            if node in seed_nodes:
                node_labels[node] = str(node)
            elif node in highlight_nodes:
                node_labels[node] = str(node)
            else:
                node_labels[node] = ''
    else:
        node_labels = {n:str(n) for n in nodes}

    nx.set_node_attributes(G, name = 'nodeLabel', values = node_labels)

    # set color of each node ####################################################################
    # G = nx.Graph(G)
    node_to_cluster = community.best_partition(G)
    nx.set_node_attributes(G, name = 'cluster', values = node_to_cluster)
    node_to_color = {}
    node_to_color = assign_colors_to_clusters(node_to_cluster, cluster_size_cut_off = cluster_size_cut_off, color_list = None)
	
    # set the position of each node
    cluster_pos = nx.spring_layout(G)
    rad_pos = bias_position_by_partition(cluster_pos, node_to_cluster, r=0.5)
	
    xpos,ypos=zip(*rad_pos.values())
    nx.set_node_attributes(G, name = 'xpos', values = dict(zip(rad_pos.keys(),[x*1000 for x in xpos])))
    nx.set_node_attributes(G, name = 'ypos', values = dict(zip(rad_pos.keys(),[y*1000 for y in ypos])))
    
    # set the title of each node
    node_titles = [str(node[0]) + '<br/>cluster = ' + str(node[1]['cluster'])
                   for node in G.nodes(data=True)]
    node_titles = dict(zip(nodes,node_titles))
    nx.set_node_attributes(G, name = 'nodeTitle', values = node_titles)

    # create the nodes_dict with all relevant fields
    nodes_dict = [{'id':str(n),
 #                  'border_width':border_width[n],
                   'degree':G.degree(n),
                   'color':node_to_color[n],
                   'node_label':node_labels[n],
                   'node_size':node_size,
                   'node_shape':node_to_shape[n],
                   'title':node_titles[n],
                   'x':rad_pos[n][0]*node_spacing,
                   'y':rad_pos[n][1]*node_spacing} for n in nodes]

    # map nodes to indices for source/target in edges
    node_map = dict(zip(nodes,range(len(nodes))))

    # create the edges_dict with all relevant fields
    edges_dict = [{'source':node_map[edges[i][0]],
                   'target':node_map[edges[i][1]],
                   'color':'grey'} for i in range(len(edges))]

    # set node_size_multiplier to increase node size as graph gets smaller
    if 'node_size_multiplier' not in kwargs.keys():
        if len(nodes) > 500:
            kwargs['node_size_multiplier'] = 3
        elif len(nodes) > 200:
            kwargs['node_size_multiplier'] = 5
        else:
            kwargs['node_size_multiplier'] = 7

    kwargs['physics_enabled'] = physics_enabled

    # if node hovering color not set, set default to black
    if 'node_color_hover_background' not in kwargs.keys():
        kwargs['node_color_hover_background'] = 'black'

    # node size determined by size in nodes_dict, not by id
    if 'node_size_field' not in kwargs.keys():
        kwargs['node_size_field'] = 'node_size'

    # node label determined by value in nodes_dict
    if 'node_label_field' not in kwargs.keys():
        kwargs['node_label_field'] = 'node_label'

    # export the network to JSON for Cytoscape
    if export_network:
        node_colors = map_node_to_color(G,'node_heat',True)
        nx.set_node_attributes(G, name = 'nodeColor', values = node_colors)
        edge_colors = map_edge_to_color(G,'edge_weight',True)
        nx.set_edge_attributes(G, name = 'edgeColor', values = edge_colors)
        visJS_module.export_to_cytoscape(G = G,export_file = export_file)

    return visJS_module.visjs_network(nodes_dict,edges_dict,**kwargs), cluster_pos, rad_pos, node_to_cluster