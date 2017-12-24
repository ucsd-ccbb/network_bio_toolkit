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


def tr_pvalues(DG_TF, DG_universe, DEG_list):

    """
        Our p-value function calculates the log of the p-value for every TF in the graph using [scipy.stats.hypergeom.logsf]
        (https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.hypergeom.html). These values help us
        determine which TF's are actually associated with our DEG's. If a TF is given a high value (because we are
        working with logs, not straight p-values), then it is likely that there is correlation between that TF and its
        DEG targets. Therefore, it is likely that TF is responsible for some of our observed gene expression.
        Note that if a TF is given a value of zero, that means none of the TF's targets were DEG's.

        Args:
            DG_TF: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes
            DG_universe: a networkx graph containing all interactions in our universe
            DEG_list: list of strings, your list of differentially expressed genes

        Returns: A sorted Pandas Series that maps a transcription factor's gene symbol to its calculated p-vlaue log.

    """

    source_nodes = list(set(zip(*DG_TF.edges())[0]))  # identifying unique source nodes in graph
    DG_universe_edges = list(DG_universe.edges())
    background_list = list(set(zip(*DG_universe_edges)[0]) | set(zip(*DG_universe_edges)[1]))

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

        TR_to_pvalue = pd.Series(TR_to_pvalue).sort_values(ascending=False)

    return TR_to_pvalue


# MJW to force use of unbiased calculation, set bias filter to 1
# MJW to force use of biased calculation, set bias filter to -1
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
        return bias_corrected_tf_zscore(DG, DEG_list, bias)
    else:
        return not_bias_corrected_tf_zscore(DG, DEG_list)


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


def rank_and_score_df(series, genes_to_rank, value_name = 'z-score', abs_value = True, act = False, remove_dups = False):

    """
        This function ranks an input set of genes (genes_to_rank) based on each gene's zscore relative to all calculated
        z-scores(series). By default, this function assumes it is sorting by z-score, but this funciton can also be
        used for pvalue. You can choose to sort by most activating (abs_value = False, act = True), most inhibiting
        (abs_value = False, act = False), or highest absolute zscore (abs_value = True).

        Args:
            series: a Pandas Series that maps all TF's to their calculated z-scores (or p-values)
            genes_to_rank: a list of strings, the genes symbols that you wish to know the rankings of
            value_name: string, this will be the title of the third column of the output Dataframe (should match what
                series maps to)
            abs_value: Boolean, True to rank based on absolute z-score, False otherwise
            act: Boolean, True to rank by most positive z-score, False to rank by most negative. Ignored if abs_value
                is True.
            remove_dups: Boolean, True to give genes with the same z-score the same rank, False to assign every gene a
                different rank even if they have the same z-score.

        Returns: A sorted Pandas Dataframe displaying each gene symbol alongside its rank and z-score

    """

    genes_to_rank = list(set(genes_to_rank)) # cannot have any duplicates
    scores = series.loc[genes_to_rank]
    ranks, num_ranks = rank(series, genes_to_rank, abs_value, act, remove_dups)

    df = pd.concat([ranks, scores], axis=1)
    df = df.rename(columns={0: 'rank', 1: value_name})

    df = df.sort_values(['rank'], ascending=True)
    return df




def top_values(series, act = True, abs_value = False, top = 10):

    """
        This function returns a sorted Pandas Series of the top (number indicated by top) genes based off z-score
        (given by series).

        Args:
            series: a Pandas Series that maps all TF's to their calculated z-scores (or p-values)
            act: Boolean, True to sort by most positive z-score, False to sort by most negative. Ignored if abs_value
                is True
            abs_value: Boolean, True to sort by absolute z-score, False otherwise
            top: the number of genes you wish to have returned

        Returns: A sorted Pandas Dataframe of the top genes, mapping each gene to its z-score

    """

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

    """
        This function ranks an input set of genes (genes_to_rank) based on each gene's z-score relative to all calculated
        z-scores(series). By default, this function assumes it is sorting by z-score, but this funciton can also be
        used for pvalue. You can choose to sort by most activating (abs_value = False, act = True), most inhibiting
        (abs_value = False, act = False), or highest absolute zscore (abs_value = True).

        Args:
            series: a Pandas Series that maps all TF's to their calculated z-scores (or p-values)
            genes_to_rank: a list of strings, the genes symbols that you wish to know the rankings of
            abs_value: Boolean, True to rank based on absolute z-score, False otherwise
            act: Boolean, True to rank by most positive z-score, False to rank by most negative. Ignored if abs_value
                is True.
            remove_dups: Boolean, True to give genes with the same z-score the same rank, False to assign every gene a
                different rank even if they have the same z-score.

        Returns: A sorted Pandas Series mapping each gene symbol to its rank

    """

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