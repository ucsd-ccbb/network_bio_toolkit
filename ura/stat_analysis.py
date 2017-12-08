"""
-------------------------------------------
Author: Mikayla Webster (m1webste@ucsd.edu)
Date: 10/13/17
-------------------------------------------
"""

import scipy
from scipy import stats
import math
import pandas as pd


def tr_pvalues(DG, db_edges, DEG_list):

    """
        Our p-value function calculates the log of the p-value for every TF in the graph using [scipy.stats.hypergeom.logsf]
        (https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.hypergeom.html). These values help us
        determine which TF's are actually associated with our DEG's. If a TF is given a high value (because we are
        working with logs, not straight p-values), then it is likely that there is correlation between that TF and its
        DEG targets. Therefore, it is likely that TF is responsible for some of our observed gene expression.
        Note that if a TF is given a value of zero, that means none of the TF's targets were DEG's.

        Args:
            DG: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes
            db_edges: list of strings, list of all genes in your experiment's universe
            DEG_list: list of strings, your list of differentially expressed genes

        Returns: A dictionary that maps a transcription factor's gene symbol to its calculated p-vlaue log.

    """

    source_nodes = list(set(zip(*DG.edges())[0]))  # identifying unique source nodes in graph
    background_list = list(set(zip(*db_edges)[0]) | set(zip(*db_edges)[1]))

    TR_to_pvalue = {}
    for TR in source_nodes:
        x = len(list(set(DG.neighbors(TR)) & set(DEG_list)))  # per TR, observed overlap between TR neighbors and DEG_list
        M = len(background_list)  # num unique nodes in universe, aka background network (STRING)
        n = len(DG.neighbors(TR))  # per TR, number of targets for that TR
        N = len(list(set(background_list) & set(DEG_list)))  # number of DEG, picked from universe "at random"

        if x == 0:
            TR_to_pvalue[TR] = 0
        else:
            TR_to_pvalue[TR] = -(scipy.stats.hypergeom.logsf(x, M, n, N, loc=0))  # remove unnecessary negative sign
            
            
        # ---------------------------------------------------------------
        # SBR: We are getting a lot of infs --> look into why this is happening
        # ---------------------------------------------------------------
        TR_to_pvalue = pd.Series(TR_to_pvalue).sort_values(ascending=False) # SBR added sorting to output

    return TR_to_pvalue

# to force use of unbiased calculation, set bias filter to 1
# to force use of biased calculation, set bias filter to -1
def tr_zscore(DG, DEG_list, bias_filter = 0.25):

    bias = calculate_bias(DG)
    if abs(bias) > bias_filter:
        print 'Graph has bias of ' + str(bias) + '. Adjusting z-score calculation accordingly.'
        return bias_corrected_tr_zscore(DG, DEG_list, bias)
    else:
        return not_bias_corrected_tr_zscore(DG, DEG_list)


def not_bias_corrected_tr_zscore(DG, DEG_list):

    """
        The goal of our z-score function is to predict the activation states of the TF's. We observe how a TF relates
        to each of its targets to make our prediction. We compare each targets' observed gene regulation (either up or
        down) and each TF-target interaction (whether it is activating or inhibiting) to conclude whether a TF is
        activating or inhibiting. A positive value indicates activating while a negative value indicates inhibiting.
        A value of zero means that we did not have enough information about the target or TF-target interaction to
        make the prediction.

        Args:
            DG: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes
            DEG_list: list of strings, your list of differentially expressed genes

        Returns: A dictionary that maps a transcription factor's gene symbol to its calculated z-score.

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
                    up_down_of_n = (DG.node[n]['updown'] / abs(DG.node[n]['updown']))

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

# ---------------------------------------------------------------
# SBR: I might suggest consolidating bias_corrected and not_bias_corrected functions
# ---------------------------------------------------------------
def bias_corrected_tr_zscore(DG, DEG_list, bias):
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


def rank_and_score_df(series, genes_to_rank, value_name = 'z-score', abs_value = True, act = False, remove_dups = False):

    scores = series.loc[genes_to_rank]
    ranks, num_ranks = rank(series, genes_to_rank, abs_value, act, remove_dups)

    df = pd.concat([ranks, scores], axis=1)
    df = df.rename(columns={0: 'rank', 1: value_name})

    df = df.sort_values(['rank'], ascending=True)
    return df




def top_values(series, act = True, abs_value = False, top = 10):

    # top activating and inhibiting, sort by strongest zscore or log(pvalue)
    if abs_value == True:
        top_series_abs = series.abs().sort_values(ascending=False).head(top)
        top_genes = list(top_series_abs.index)
        top_values = [series[gene] for gene in top_genes]
        return pd.Series(top_values, index=top_genes)

    # top activating
    if act == True:
        return series.sort_values(ascending=False).head(top)

    # top inhibiting
    else:
        return series.sort_values(ascending=True).head(top)


def rank(series, genes_to_rank, abs_value=False, act=True, remove_dups = False):

    # genes with the same p_value/zscore will be given different ranks depending on the order they are stored in
    if remove_dups == False:
        if abs_value == True:
            sorted_dict = sorted(dict(series.abs()).items(), key=lambda x: x[1], reverse=True)
        else:
            if act == True:
                sorted_dict = sorted(dict(series).items(), key=lambda x: x[1], reverse=True)
            else:
                sorted_dict = sorted(dict(series).items(), key=lambda x: x[1], reverse=False)

        genes = zip(*sorted_dict)[0]
        index = range(len(sorted_dict))
        gene_to_index = dict(zip(genes, index)) # mapping genes to their order in the sorted gene list

        return_series = pd.Series({k: gene_to_index.get(k, None) for k in genes_to_rank}).sort_values(ascending=True)
        num_ranks = len(gene_to_index)
        return return_series, num_ranks

    # genes with the same p_value will be given the same rank
    else:
        if abs_value == True:
            sorted_dict = sorted(dict(series.abs()).items(), key=lambda x: x[1], reverse=True)
            rank_values = sorted(set(abs(series.values)), reverse = True)
        else:
            if act == True:
                sorted_dict = sorted(dict(series).items(), key=lambda x: x[1], reverse=True)
                rank_values = sorted(set(series.values), reverse = True)
            else:
                sorted_dict = sorted(dict(series).items(), key=lambda x: x[1], reverse=False)
                rank_values = sorted(set(series.values), reverse = False)

        genes = zip(*sorted_dict)[0]
        index = range(len(rank_values))
        value_to_index = dict(zip(rank_values, index)) # mapping the sorted set of possible values to their order

        if abs_value == True:
            gene_index_list = [value_to_index[abs(series[gene])] for gene in genes]  # finding the rank for each gene
        else:
            gene_index_list = [value_to_index[series[gene]] for gene in genes] # finding the rank for each gene

        gene_to_index = dict(zip(genes, gene_index_list)) # mapping genes to their rank

        return_series = pd.Series({k: gene_to_index.get(k, None) for k in genes_to_rank}).sort_values(ascending=True)
        num_ranks = str(len(index))

        return return_series, num_ranks