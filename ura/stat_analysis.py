"""
-------------------------------------------
Author: Mikayla Webster (m1webste@ucsd.edu)
Date: 10/13/17
-------------------------------------------
"""

import scipy
from scipy import stats
import math


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
            DEG_list: list of strings, your list of differencially expressed genes

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
            TR_to_pvalue[TR] = 0;
        else:
            TR_to_pvalue[TR] = -(scipy.stats.hypergeom.logsf(x, M, n, N, loc=0))  # remove unnecessary negative sign

    return TR_to_pvalue


def tr_zscore(DG, DEG_list):

    """
        The goal of our z-score function is to predict the activation states of the TF's. We observe how a TF relates
        to each of its targets to make our prediction. We compare each targets' observed gene regulation (either up or
        down) and each TF-target interaction (whether it is activating or inhibiting) to conclude whether a TF is
        activating or inhibiting. A positive value indicates activating while a negative value indicates inhibiting.
        A value of zero means that we did not have enough information about the target or TF-target interaction to
        make the prediction.

        Args:
            DG: Digraph, a directed networkx graph with edges mapping from transcription factors to expressed genes
            DEG_list: list of strings, your list of differencially expressed genes

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

    return TR_to_zscore


def calculate_bias(DG, DEG_list):
    source_nodes = list(set(zip(*DG.edges())[0]))  # identifying unique source nodes in graph

    TR_to_bias = {}
    for TR in source_nodes:

        N_up = 0  # number of up regulated target
        N_down = 0  # number of down regulated targets

        N_act = 0  # number of activating edges
        N_inh = 0  # number of inhibiting edges

        N_problem = 0  # number of edges with errorous calculations

        TRs_DEG_neighbors = set(DG.neighbors(TR)) & set(DEG_list)

        for n in TRs_DEG_neighbors:

            if ((str(type(DG)) == '<class \'networkx.classes.multidigraph.MultiDiGraph\'>') | (
                str(type(DG)) == '<class \'networkx.classes.multigraph.MultiGraph\'>')):
                for i in range(len(DG[TR][n])): # have to take into account multiple edges with the same mapping
                    # count up edge signs
                    sign_of_edge = DG[TR][n][i]['sign']
                    if sign_of_edge == 1:
                        N_act += 1
                    elif sign_of_edge == -1:
                        N_inh += 1
                    else:
                        N_problem += 1
                        print "Issue with edge (" + str(TR) + ',' + str(n) + ') A/I'

                    # count up node regulations
                    up_down_of_n = (DG.node[n]['updown'] / abs(DG.node[n]['updown']))
                    if up_down_of_n == 1:
                        N_up += 1
                    elif up_down_of_n == -1:
                        N_down += 1
                    else:
                        N_problem += 1
                        print "Issue with edge (" + str(TR) + ',' + str(n) + ') up/down'

            else:
                # count up edge signs
                sign_of_edge = DG[TR][n]['sign']
                if sign_of_edge == 1:
                    N_act += 1
                elif sign_of_edge == -1:
                    N_inh += 1
                else:
                    N_problem += 1
                    print "Issue with edge (" + str(TR) + ',' + str(n) + ') A/I'

                # count up node regulations
                up_down_of_n = (DG.node[n]['updown'] / abs(DG.node[n]['updown']))
                if up_down_of_n == 1:
                    N_up += 1
                elif up_down_of_n == -1:
                    N_down += 1
                else:
                    N_problem += 1
                    print "Issue with edge (" + str(TR) + ',' + str(n) + ') up/down'

        # calculate up down bias
        if (N_up + N_down) != 0:
            u_data = (N_up - N_down) / float(N_up + N_down)
        else:
            u_data = 0

        # calculate act-inh bias
        if (N_act + N_inh) != 0:
            u_TR = (N_act - N_inh) / float(N_act + N_inh)
        else:
            u_TR = 0

        # calculate overall bias
        u = u_data * u_TR
        TR_to_bias[TR] = u

    return TR_to_bias


def bias_corrected_tr_zscore(DG, DEG_list, TR_to_bias):
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

        u = TR_to_bias[TR]

        # calculate bias-corrected z-score
        z_score_top = 0
        z_score_bottom = 0
        for i in range(len(w)):
            z_score_top += (w[i] * (x[i] - u))
            z_score_bottom += (w[i] * w[i])
        z_score = z_score_top / ((z_score_bottom) ** (1 / 2))

        TR_to_zscore[TR] = z_score

    return TR_to_zscore