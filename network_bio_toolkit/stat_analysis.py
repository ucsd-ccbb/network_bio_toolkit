"""
-------------------------------------------
Author: Mikayla Webster (13webstermj@gmail.com)
Date: 10/13/17
-------------------------------------------
"""

import scipy
import math
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import seaborn as sns

import visJS2jupyter.visJS_module as visJS_module # "pip install visJS2jupyter"
import create_graph # from URA package


# --------------------- P-VALUE FUNCTIONS ---------------------------#


def tf_pvalues(DG_TF, DG_universe, DEG_list):

    """
        Our p-value function calculates the log of the p-value for every TF in the graph using [scipy.stats.hypergeom.logsf]
        (https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.hypergeom.html). These values help us
        determine which TF's are actually associated with our DEG's. If a TF is given a high value (because we are
        working with logs, not straight p-values), then it is likely that there is correlation between that TF and its
        DEG targets. Therefore, it is likely that TF is responsible for some of our observed gene expression.
        Note that if a TF is given a value of zero, that means none of the TF's targets were DEG's.

        Args:
            DG_TF: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes (filtered)
            DG_universe: a networkx graph containing all interactions in our universe (not filtered)
            DEG_list: list of strings, your list of differentially expressed genes

        Returns: A sorted Pandas Series that maps a transcription factor's gene symbol to its calculated p-vlaue log.

    """

    source_nodes = list(set(zip(*DG_TF.edges())[0]))  # identifying unique source nodes in graph
    background_list = list(DG_universe.nodes()) # list of all unique nodes in universe

    TR_to_pvalue = {}
    x_n_to_p_score = {}
    M = len(background_list)  # num unique nodes in universe, aka background network (STRING)
    N = len(list(set(background_list) & set(DEG_list)))  # number of DEG, picked from universe "at random"

    for TR in source_nodes:
        x = len(list(set(DG_TF.neighbors(TR)) & set(DEG_list)))  # per TR, observed overlap between TR neighbors and DEG_list
        n = len(DG_TF.neighbors(TR))  # per TR, number of targets for that TR

        if (x,n) in x_n_to_p_score: # if we have computed this value before
            TR_to_pvalue[TR] = x_n_to_p_score[(x,n)]

        else:
            if x == 0:
                TR_to_pvalue[TR] = 0
            elif x == n:
                TR_to_pvalue[TR] = float('Inf')
            else:
                TR_to_pvalue[TR] = -(scipy.stats.hypergeom.logsf(x, M, n, N, loc=0))  # remove unnecessary negative sign

            x_n_to_p_score[(x,n)] = TR_to_pvalue[TR] # record that we have calculated this value

        TR_to_pvalue = pd.Series(TR_to_pvalue).sort_values(ascending = False)

    return TR_to_pvalue



def tf_target_enrichment(DG_TF, DG_universe, DEG_list):

    """
        Sets up tf_pvalues() to calculate enrichment of transcription factors and their
        targets.
    """

    to_return = tf_pvalues(DG_TF, DG_universe, DEG_list)
    to_return = to_return.rename('tf-target enrichment')
    return to_return


def tf_enrichment(TF_list, DEG_full_graph, DEG_list):

    """
        Sets up tf_pvalues() to calculate enrichment of transcription factors.
    """

    # make a graph that will manipulate our p-value function to calculate one TF-enrichment p-value
    dummy = 'TF_ENRICHMENT'
    dummy_list = [dummy] * len(TF_list)
    edges = zip(dummy_list, TF_list)
    G = nx.DiGraph(edges)

    return tf_pvalues(G, DEG_full_graph, DEG_list)



# --------------------- Z-SCORE FUNCTIONS ---------------------------#

	
def tf_zscore(DG, DEG_list, bias_filter = 0.25):

    """
        The goal of our z-score function is to predict the activation states of the TF's. We observe how a TF relates
        to each of its targets to make our prediction. We compare each targets' observed gene regulation (either up or
        down) and each TF-target interaction (whether it is activating or inhibiting) to conclude whether a TF is
        activating or inhibiting. A positive value indicates activating while a negative value indicates inhibiting.
        A value of zero means that we did not have enough information about the target or TF-target interaction to
        make the prediction.

        This function call one of two helper z-score functions, either bias_corrected_tf_zscore or not_bias_corrected_tf_zscore,
        based on how biased the graph is (indicated by the bias_filter parameter). The "bias" of the graph is a number that
        indicates if the graph has notibly more activating or inhibiting edges, and to what degree. It is calculated using our
        calculate_bias function.

        **If the user wishes to explicitly use the biased z-score formula (bias_corrected_tf_zscore), set bias_filter to 0.
        For the unbiased formula (not_bias_corrected_tf_zscore), set bias_filter to 1.

        Args:
            DG: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes
                ** DG must have 'updown' and 'sign' attributes **
            DEG_list: list of strings, your list of differentially expressed genes
            bias_filter: number between 0 and 1, threshold to calculate z-score using biased formula

        Returns: A sorted Pandas Series that maps a transcription factor's gene symbol to its calculated z-score.

    """

    bias = calculate_bias(DG)
    if abs(bias) > bias_filter:
        print 'Graph has bias of ' + str(bias) + '. Adjusting z-score calculation accordingly.'
        to_return = bias_corrected_tf_zscore(DG, DEG_list, bias)
        to_return = to_return.rename('biased-calculated z-scores')
        return to_return
    else:
        to_return = not_bias_corrected_tf_zscore(DG, DEG_list)
        to_return = to_return.rename('unbiased-calculated z-scores')
        return to_return


def not_bias_corrected_tf_zscore(DG, DEG_list):

    """
        The goal of our z-score function is to predict the activation states of the TF's. We observe how a TF relates
        to each of its targets to make our prediction. We compare each targets' observed gene regulation (either up or
        down) and each TF-target interaction (whether it is activating or inhibiting) to conclude whether a TF is
        activating or inhibiting. A positive value indicates activating while a negative value indicates inhibiting.
        A value of zero means that we did not have enough information about the target or TF-target interaction to
        make the prediction.

        This calculation does NOT take into account network bias.

        Args:
            DG: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes.
                ** DG must have 'updown' and 'sign' attributes **
            DEG_list: list of strings, your list of differentially expressed genes

        Returns: A sorted Pandas Series that maps a transcription factor's gene symbol to its calculated z-score.

    """

    source_nodes = list(set(zip(*DG.edges())[0]))  # identifying unique source nodes in graph

    TR_to_zscore = {}
    for TR in source_nodes:
        N_minus = 0  # number of inhibiting predicting DEG edges
        N_plus = 0  # number of activating predicting DEG edges
        N_zero = 0  # number of edges with errorous calculations

        TRs_DEG_neighbors = set(DG.neighbors(TR)) & set(DEG_list)
        for n in TRs_DEG_neighbors:

            if ((str(type(DG)) == '<class \'networkx.classes.multidigraph.MultiDiGraph\'>') | (str(type(DG)) == '<class \'networkx.classes.multigraph.MultiGraph\'>')):
                for i in range(len(DG[TR][n])): # have to take into account multiple edges with the same mapping
                    sign_of_edge = DG[TR][n][i]['sign']
					
                    if float(DG.node[n]['updown']) != 0:
                        up_down_of_n = (DG.node[n]['updown'] / float(abs(DG.node[n]['updown'])))
                    else:
                        print "Edge (" + str(TR) + ',' + str(n) + ') has updown of zero'

                    # predict whether this neighbor thinks the TR is Act. or Inhib.
                    if ((sign_of_edge * up_down_of_n) == 1):
                        N_plus += 1
                    elif ((sign_of_edge * up_down_of_n) == -1):
                        N_minus += 1
                    else:
                        N_zero += 1  # mark an error if could not predict
                        print "Issue with edge (" + str(TR) + ',' + str(n) + ')'
            else:
                sign_of_edge = DG[TR][n]['sign']
                up_down_of_n = (DG.node[n]['updown'] / abs(DG.node[n]['updown']))

                # predict whether this neighbor thinks the TR is Act. or Inhib.
                if ((sign_of_edge * up_down_of_n) == 1):
                    N_plus += 1
                elif ((sign_of_edge * up_down_of_n) == -1):
                    N_minus += 1
                else:
                    N_zero += 1  # mark an error if could not predict
                    print "Issue with edge (" + str(TR) + ',' + str(n) + ')'

        if N_zero != 0:
            print "Could not attribute activated or inhibiting trait to " + str(N_zero) + 'nodes'

        # prevent a divide-by-zero calculation
        N = N_plus + N_minus
        if N == 0:
            z_score = 0
        else:
            z_score = (N_plus - N_minus) / float(math.sqrt(N))

        TR_to_zscore[TR] = z_score
        # create zscore dict where positive means activating
        # negative means inhibiting
        # 0 means could not be calculated

        TR_to_zscore = pd.Series(TR_to_zscore).sort_values(ascending=False) # SBR sorted output
        
    return pd.Series(TR_to_zscore)


def calculate_bias(DG):

    """
        The goal of our z-score function is to predict the activation states of the TF's. We observe how a TF relates
        to each of its targets to make our prediction. We compare each targets' observed gene regulation (either up or
        down) and each TF-target interaction (whether it is activating or inhibiting) to conclude whether a TF is
        activating or inhibiting. A positive value indicates activating while a negative value indicates inhibiting.
        A value of zero means that we did not have enough information about the target or TF-target interaction to
        make the prediction.

        This calculation does NOT take into account network bias.

        Args:
            DG: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes.
                ** DG must have 'updown' and 'sign' attributes **

        Returns: A number between 0 and 1 that is the graph's bias

    """

    # calculate bias for up downs
    data = list(zip(*DG.nodes(data=True))[1])
    ups = [dict_list['updown'] for dict_list in data if dict_list['updown'] > 0]
    downs = [dict_list['updown'] for dict_list in data if dict_list['updown'] < 0]
    N_up = len(ups)
    N_down = len(downs)

    if (N_up + N_down) != 0:
        u_data = (N_up - N_down) / float(N_up + N_down)
    else:
        u_data = 0

    # calculate bias for activating and inhbiting
    data = list(zip(*DG.edges(data=True))[2])
    act = [dict_list['sign'] for dict_list in data if dict_list['sign'] == 1]
    inh = [dict_list['sign'] for dict_list in data if dict_list['sign'] == -1]
    N_act = len(act)
    N_inh = len(inh)

    if (N_act + N_inh) != 0:
        u_TR = (N_act - N_inh) / float(N_act + N_inh)
    else:
        u_TR = 0

    return u_data * u_TR

def bias_corrected_tf_zscore(DG, DEG_list, bias):

    """
        The goal of our z-score function is to predict the activation states of the TF's. We observe how a TF relates
        to each of its targets to make our prediction. We compare each targets' observed gene regulation (either up or
        down) and each TF-target interaction (whether it is activating or inhibiting) to conclude whether a TF is
        activating or inhibiting. A positive value indicates activating while a negative value indicates inhibiting.
        A value of zero means that we did not have enough information about the target or TF-target interaction to
        make the prediction.

        This calculation assumes background network bias, indicated by argument bias.

        Args:
            DG: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes
                ** DG must have 'updown' and 'sign' attributes **
            DEG_list: list of strings, your list of differentially expressed genes
            bias: a number between 0 and 1 that is the background network's bias

        Returns: A sorted Pandas Series that maps a transcription factor's gene symbol to its calculated z-score.

    """

    source_nodes = list(set(zip(*DG.edges())[0]))  # identifying unique source nodes in graph

    TR_to_zscore = {}
    for TR in source_nodes:
        w = []
        x = []
        TRs_DEG_neighbors = set(DG.neighbors(TR)) & set(DEG_list)
        for n in TRs_DEG_neighbors:

            if ((str(type(DG)) == '<class \'networkx.classes.multidigraph.MultiDiGraph\'>') | (
                        str(type(DG)) == '<class \'networkx.classes.multigraph.MultiGraph\'>')):
                for i in range(len(DG[TR][n])):  # have to take into account multiple edges with the same mapping
                    sign_of_edge = DG[TR][n][i]['sign']
                    up_down_of_n = (DG.node[n]['updown'] / abs(DG.node[n]['updown']))

                    # predict whether this neighbor thinks the TR is Act. or Inhib.
                    prediction = sign_of_edge * up_down_of_n
                    if ((prediction == 1) | (prediction == -1)):
                        x.append(prediction)
                    else:
                        print "Issue with edge (" + str(TR) + ',' + str(n) + ')'

                        # keep track of each target's weight
                    w.append(DG[TR][n][i]['weight'])

            else:
                sign_of_edge = DG[TR][n]['sign']
                up_down_of_n = (DG.node[n]['updown'] / abs(DG.node[n]['updown']))

                # predict whether this neighbor thinks the TR is Act. or Inhib.
                prediction = sign_of_edge * up_down_of_n
                if ((prediction == 1) | (prediction == -1)):
                    x.append(prediction)
                else:
                    print "Issue with edge (" + str(TR) + ',' + str(n) + ')'

                    # keep track of each target's weight
                w.append(DG[TR][n]['weight'])

        # calculate bias-corrected z-score
        z_score_top = 0
        z_score_bottom = 0
        for i in range(len(w)):
            z_score_top += (w[i] * (x[i] - bias))
            z_score_bottom += (w[i] * w[i])
        z_score = z_score_top / ((z_score_bottom) ** (1 / 2))

        TR_to_zscore[TR] = z_score
        
        TR_to_zscore = pd.Series(TR_to_zscore).sort_values(ascending=False) # SBR sorted output

    return pd.Series(TR_to_zscore)


# --------------------- DISPLAY FUNCTIONS ---------------------------#



def top_values(z_score_series, DEG_to_pvalue, DEG_to_updown, act = True, abs_value = False, top = 10):

    """
        This function returns a sorted Pandas Series of the top (number indicated by top) genes based off z-score
        (given by z_score_series).

        Args:
            z_score_series: a Pandas Series that maps all TF's to their calculated z-scores (or p-values)
            DEG_to_pvalue: Dict, maps gene names to their (adj) p-values found in the expression file
            DEG_to_updown: Dict, maps gene names to their (log) fold change found in the expression file
            act: Boolean, True to sort by most positive z-score, False to sort by most negative. Ignored if abs_value
                is True
            abs_value: Boolean, True to sort by absolute z-score, False otherwise
            top: the number of genes you wish to have returned

        Returns: A sorted Pandas Dataframe of the top genes, mapping each gene to its z-score
    """
    
    # top activating and inhibiting, sort by strongest zscore or log(pvalue)
    if abs_value == True:
        top_series_abs = z_score_series.abs().sort_values(ascending=False).head(top)
        top_genes = list(top_series_abs.index)
        top_values = [z_score_series[gene] for gene in top_genes]
        z_score_series = pd.Series(top_values, index=top_genes, name ='z-score')

    # top activating
    if act == True:
        z_score_series = z_score_series.sort_values(ascending=False).head(top)
        z_score_series = z_score_series.rename('z-score')

    # top inhibiting
    else:
        z_score_series = z_score_series.sort_values(ascending=True).head(top)
        z_score_series = z_score_series.rename('z-score')

    # list of our top genes
    top_genes = list(z_score_series.index)

    # get p-value series
    TF_to_adjp = {gene:DEG_to_pvalue[gene] if gene in DEG_to_pvalue else None for gene in top_genes}
    p_series = pd.Series(TF_to_adjp, index=top_genes, name = '(adj) p-value')

    # get fld change series
    TF_to_foldchange = {gene: DEG_to_updown[gene] if gene in DEG_to_updown else None for gene in top_genes}
    fld_series = pd.Series(TF_to_foldchange, index=top_genes, name = '(log) fold change')

    return pd.concat([z_score_series, p_series, fld_series], axis=1)




def compare_genes(z_scores, genes_to_rank, fig_size = (12, 7), font_size = 12, anno_vert_dist = 0.025):

    """
        Plots a histogram of the distribution of z-scores, and marks the genes indicated in genes_to_rank
        so that their values can be visually compared.

        Args:
            z_scores: List, z-scores of all genes in your analysis
            genes_to_rank: List, genes you wish to highlight in the distribution
            fig_size: (Int, Int), the dimensions of the histogram.
            font_size: Int, the font size of the genes_to_rank labels
            anno_vert_dist: Float, the distance separating the genes_to_rank labels from the x-axis

        Returns: A seaborn based histogram.
    """

    # We don't want to plot the zero z-scores
    z_scores_hist = [x for x in z_scores if x != 0]

    # group together gene names based on their z-scores
    annotate_list = []
    new_genes_to_rank = []
    for i in range(len(genes_to_rank)):
        if genes_to_rank[i] in z_scores:
            new_genes_to_rank.append(genes_to_rank[i])
            annotate_list = put_in_bucket(z_scores, annotate_list, genes_to_rank[i])
        else:
            print str(genes_to_rank[i]) + ' is not a valid transcription factor in our graph.'
    genes_to_rank = new_genes_to_rank

    # base x-axis off z-scores
    x_points = []
    for i in range(len(annotate_list)):
        x_points.append(z_scores[annotate_list[i][0]])

    y_points = [0] * len(annotate_list)

    # plot points
    plt.figure(figsize = fig_size)
    ax = sns.distplot(z_scores_hist, kde = True)
    plt.scatter(x_points, y_points, marker = '^', s = 200, c = 'r')

    # annotate points
    for i in range(len(annotate_list)):
        to_print = str(annotate_list[i]).replace("[", "")
        to_print = to_print.replace("]", "")
        to_print = to_print.replace("\'", "")
        ax.annotate(to_print,
                    xy = (x_points[i],
                        anno_vert_dist),
                    rotation = 90,
                    horizontalalignment = 'center',
                    verticalalignment = 'bottom',
                    fontsize=font_size)

		
def vis_tf_network(DG, tf, DEG_list,
                   DEG_to_pvalue, DEG_to_updown,
                   directed_edges = False,
                   node_spacing = 2200,
                   color_non_DEGs = False,
                   color_map = plt.cm.bwr,
				   graph_id = 0,
                   tf_size_amplifier = 8,
				   alpha = 1.0, 
				   color_vals_transform = None,
				   ceil_val=10,
                   color_max_frac = 1.0,
				   color_min_frac = 0.0,
				   vmin = None,
				   vmax = None,
				   tf_shape = 'star'
                   ):
				   
    """
        This fuction visualizes the network consisting of one transcription factor and its downstream regulated genes. Nodes are colored with 
        a shade of blue or red calculated based on the fold change of each gene given in the file DEG_filename. Red is up-regulated. Blue is 
        down-regulated. White nodes did not have enough information. Darker red/blue indicates a stronger (larger absolute value) fold change 
        value. Node size is deterined by adjusted p-value, also from DEG_filename. Larger nodes have a more significant p-value. Nodes from 
        DEG_list will have a black border. Activating edges are red. Inhibiting edges are blue. 

        Args:
            DG: Networkx graph, your (potentially TF-filtered) background networkx
            tf: String, all caps gene symbol of the regulator whose sub-network you wish to visualizes
            DEG_list: List, list of DEG's as strings, output of create_graph.creat_DEG_list
            DEG_to_pvalue: Dict, maps gene names to their (adj) p-values found in the expression file
            DEG_to_updown: Dict, maps gene names to their (log) fold change found in the expression file
            directed_edges: Boolean, True to include directional arrows on edges, False otherwise
            node_spacing: Int, increase this number if your nodes are too close together (for a graph with many nodes) or decrease 
                if they are too far apart (for a graph with fewer nodes)
            color_non_DEGs: Boolean, True will color all nodes in graph with their log fold change, False will leave non-DEG genes grey.
            color_map: matplotlib.cm.*, a matplotlib colormap
            graph_id: Int, change this number to display multiple graphs in one notebook
            tf_size_amplifier: Int, size multiplier applied only to the transcription factor
            alpha: Float, alpha parameter to visJS_module.return_node_to_color()
            color_vals_transform: String, function name. Either 'log', 'sqrt', 'ceil', or 'None'. 
                color_vals_transform parameter to visJS_module.return_node_to_color()
            ceil_val: Float, value cut-off for data. ceil_val parameter to visJS_module.return_node_to_color()
            color_max_frac: Float, color_max_frac parameter to visJS_module.return_node_to_color()
            color_min_frac: Float, color_min_frac parameter to visJS_module.return_node_to_color()
            vmin: Float, vmin parameter to visJS_module.return_node_to_color()
            vmax: Float, vmax parameter to visJS_module.return_node_to_color()
            tf_shape: String, ex. 'star', 'circle', 'dot'... Shape to set transcription factor

        Returns: HTML output that will display an interactive network in a jupyter notebooks cell.

    """
	
    DEG_list = list(DEG_list)
    
    #------------ GRAPH --------------#

    # create sub graph
    node_list = DG.neighbors(tf) + [tf]
    G = DG.subgraph(node_list) #get only top z-score node and neighbors
    
    # remove arrows pointing to tf
    if directed_edges == True:
        tf_in_edges = G.in_edges(tf)
        G.remove_edges_from(tf_in_edges)
    
    # define nodes and edges
    nodes = list(G.nodes())
    edges = list(G.edges()) 
    
    #------------ NODES --------------#

    # define the initial positions of the nodes using networkx's spring_layout function
    pos = nx.spring_layout(G)
    
    # define node colors
    node_to_fld = {n: DEG_to_updown[n] for n in nodes if n in DEG_to_updown} # keep only those in graph G
    nx.set_node_attributes(G, 'fold_change', 0) # give all nodes a default fold change of zero
    nx.set_node_attributes(G, 'fold_change', node_to_fld) # overwrite with actual fold change for the nodes that have one
    node_to_color = visJS_module.return_node_to_color(G, field_to_map='fold_change', 
                                                        cmap = color_map, 
                                                        alpha = alpha, 
                                                        color_vals_transform = color_vals_transform,
                                                        ceil_val = ceil_val,
                                                        color_max_frac = color_max_frac,
                                                        color_min_frac = color_min_frac,
                                                        vmin = vmin,
                                                        vmax = vmax)

    if color_non_DEGs == False: # if they don't want to see non-DEG's colors
        grey_list = [x for x in nodes if x not in DEG_list]
        for n in grey_list:
                node_to_color[n] = 'grey'
    
    # define node size
    nts = {n: (2 + (-np.log(DEG_to_pvalue[n]))) for n in nodes if n in DEG_to_updown}  # keep only those in graph G
    avg = np.average(list(nts.values()))
    max_size = np.max(list(nts.values()))
    node_to_size = {n:nts[n] if n in nts else avg for n in nodes}
    node_to_size[tf] = max_size + tf_size_amplifier
    
    # give DEG nodes an outline
    node_to_border_width = {n:2 if n in list(DEG_list) else 0 for n in nodes}
    nx.set_node_attributes(G, name = 'nodeOutline', values = node_to_border_width) # give all nodes a default fold change of zero

    # define node shape
    all_circles = ['circle']*len(nodes)
    node_to_shape = dict(zip(nodes, all_circles))
    node_to_shape[tf] = tf_shape
    nx.set_node_attributes(G, name = 'shape', values = node_to_shape)  # give all nodes a default fold change of zero
    
    # define node attributes
    nodes_dict = [{"id":n,
                   "color": node_to_color[n],
                   "border_width": node_to_border_width[n],
                   "size_field": node_to_size[n],
                   "node_shape": node_to_shape[n],
                   "x":pos[n][0] * node_spacing,
                   "y":pos[n][1] * node_spacing} for n in nodes]
    node_map = dict(zip(nodes,range(len(nodes))))  # map to indices for source/target in edges
    
    #------------ EDGES --------------#
    
    # define edge colors
    color = {-1:'blue', 1:'red', 0: 'grey'} # red = act/up, blue = inh/down, no info = grey
    edge_to_color = {edge[0:2]:color[np.sign(edge[2]['sign'])] for edge in list(G.edges(data = True))}

    # define edge attributes
    edges_dict = [{"source":node_map[edges[i][0]], 
                   "target":node_map[edges[i][1]], 
                   "color":edge_to_color[edges[i]]} for i in range(len(edges))]

    # display results
    return visJS_module.visjs_network(nodes_dict,edges_dict,
											node_color_highlight_border = 'black',
											node_color_hover_border = 'black',
                                            edge_arrow_to = directed_edges,
                                            node_size_multiplier = 10,
                                            node_size_transform = '',
                                            node_size_field = 'size_field',
                                            node_font_size = 35,
                                            node_border_width = 5,
                                            node_font_color = 'black',
                                            edge_color_hover = 'orange',
                                            edge_color_highlight = 'orange',
                                            edge_width = 1.5,
                                            edge_hoverWidth = 5,
                                            edge_selection_width = 5,
                                            node_border_width_selected = 0,
                                            graph_id = graph_id)




def to_csv(out_filename, z_score_series, DEG_to_pvalue, DEG_to_updown, tf_target_enrichment, DG_TF):

    """
        Outputs information gathered from URA Analysis to a single csv file.

        Args:
            out_filename: String, filepath to write to.
            z_score_series: Pandas series, z-scores obtained from calling tf_zscore()
            DEG_to_pvalue: Dict, maps gene names to their (adj) p-values found in the expression file
            DEG_to_updown: Dict, maps gene names to their (log) fold change found in the expression file
            tf_target_enrichment:Pandas series, p-values obtained from calling tf_pvalues()
            DG_TF: Networkx, String network filtered by transcription factors, output of create_graph.load_STRING_to_digraph()
    """
    
    top_overall = top_values(z_score_series, DEG_to_pvalue, DEG_to_updown, act = False, abs_value = True,
                             top = len(z_score_series))
    top_overall.index.name = 'TF'

    top_overall['TF_target_enrichment'] = tf_target_enrichment
    TF_to_neighbors = [list([str(x) for x in DG_TF.neighbors(tf)]) for tf in top_overall.index]
    top_overall['targets'] = TF_to_neighbors

    top_overall.to_csv(path_or_buf = out_filename, sep='\t')

    # Read in the file
    with open(out_filename, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace('[', '')
    filedata = filedata.replace(']', '')
    filedata = filedata.replace('\'', '')

    # Write the file out again
    with open(out_filename, 'w') as file:
        file.write(filedata)



# --------------------- HELPER FUNCTIONS ---------------------------#


def put_in_bucket(d, l, value):

    """
        If this helper function is called on every gene in our analysis, it will group all 
        genes with the same z-score into the same list. Each of those lists of same-z-score 
        genes is stored in l.

        Args:
            d: Dict, in the context of this file, a dict of genes to z-scores
            l: List of Lists, represents out buckets. Initiate with an empty list, the continually pass in
                l to build up your buckets 
            value: String, the current gene we would like to put into a bucket

        Returns: the list now containing value in its proper bucket within l.

    """

    # dummy list to prevent ['string'] -> ['s','t','r','i','n','g']
    dummy = []

    # if list is empty, init it
    if len(l) == 0:
        dummy.append(value)
        return dummy
    else:

        # else search to see if the value fits into an existing bucket
        for i in range(len(l)):

            # along the way, make sure this is a list of lists
            if type(l[i]) != list:
                dummy.append(l[i])
                l[i] = dummy
                dummy = []

            # aka find a bucket with same z-score as value's
            if d[l[i][0]] == d[value]:
                l[i].append(value)
                return l

        # if our value (gene) doesn't have a bucket to go in, make a new one at the end of the list
        dummy.append(value)
        l.append(dummy)
        return l