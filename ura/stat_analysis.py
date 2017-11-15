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
        Calculates a log of a p-value for every transcription factor in the graph. These values determine
        how statistically likely it that a specific transcription factor is randomly associated with the set
        of differentially expressed genes found to have been associated with it in practice.

        Args:
            DG: A directed networkX graph with edges mapping from transcription factors to expressed genes
            background_list: list of all genes in your experiment's universe
            DEG_list: your list of differencially expressed genes

        Returns: A dictionary that maps a transcription factor's gene symbol to its calculated p-vlaue log.

    """

    source_nodes = list(set(zip(*DG.edges())[0]))  # identifying unique source nodes in graph
    background_list = list(set(zip(*db_edges)[0]) | set(zip(*db_edges)[1]))

    TR_to_pvalue = {}
    for TR in source_nodes:
        x = len(list(
            set(DG.neighbors(TR)) & set(DEG_list)))  # per TR, observed overlap between TR neighbors and DEG_list
        M = len(background_list)  # num unique nodes in universe, aka background network (STRING)
        n = len(DG.neighbors(TR))  # per TR, number of targets for that TR
        N = len(list(set(background_list) & set(DEG_list)))  # number of DEG, picked from universe "at random"

        TR_to_pvalue[TR] = -(scipy.stats.hypergeom.logsf(x, M, n, N, loc=0))  # remove unnecessary negative sign

    return TR_to_pvalue


def tr_zscore(DG, DEG_list):
    source_nodes = list(set(zip(*DG.edges())[0]))  # identifying unique source nodes in graph

    TR_to_zscore = {}
    for TR in source_nodes:
        N_minus = 0  # number of inhibiting predicting DEG edges
        N_plus = 0  # number of activating predicting DEG edges
        N_zero = 0  # number of edges with errorous calculations

        TRs_DEG_neighbors = set(DG.neighbors(TR)) & set(DEG_list)
        for n in TRs_DEG_neighbors:
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

        TR_to_zscore[TR] = z_score  # create zscore dict where 1 means activating
        # -1 means inhibiting
        # 0 means could not be calculated

    return TR_to_zscore