"""
-------------------------------------------
Author: Mikayla Webster (m1webste@ucsd.edu)
Date: 2/5/18
-------------------------------------------
"""

import create_graph
import stat_analysis
import matplotlib.pyplot as plt
import copy

class Upstream:

    def __init__(self, gene_type = 'symbol'):

        if (gene_type != 'entrez') & (gene_type != 'symbol'):
            raise ValueError('Please specify either \"entrez\" or \"symbol\" as your initialization argument.\n')

        self.gene_type = gene_type

        # create_graph
        self.TF_list = None
        self.DG_universe = None
        self.DG_TF = None
        self.DEG_list = None
        self.DEG_filename = None
        self.DEG_to_pvalue = None # (adj) p-value of all genes
        self.DEG_to_updown = None # (log) fold change of all genes
        self.DEG_full_graph = None

        # stat_analysis
        self.tf_target_enrichment = None
        self.tf_enrichment = None
        self.z_scores = None

    def __repr__(self):
        return self.create_to_string(self)

    def __str__(self):
        return self.create_to_string(self)

    def create_to_string(self):
        return "TODO: create_to_string()"


    # ----------------------------- GETTERS AND SETTERS -------------------------- #
    # So our users don't accidentally modify an instance variable without realizing

    # to ask to see a variable
    def get(self, item):

        # map string to actual instance variable
        string_to_item = {}
        string_to_item['gene_type'] = self.gene_type
        string_to_item['TF_list'] = self.TF_list
        string_to_item['DG_universe'] = self.DG_universe
        string_to_item['DG_TF'] = self.DG_TF
        string_to_item['DEG_list'] = self.DEG_list
        string_to_item['DEG_filename'] = self.DEG_filename
        string_to_item['DEG_to_pvalue'] = self.DEG_to_pvalue
        string_to_item['DEG_to_updown'] = self.DEG_to_updown
        string_to_item['DEG_full_graph'] = self.DEG_full_graph
        string_to_item['tf_target_enrichment'] = self.tf_target_enrichment
        string_to_item['tf_enrichment'] = self.tf_enrichment
        string_to_item['z_scores'] = self.z_scores

        # check if the user input is valid
        try:
            to_return = string_to_item[item]
        except:
            print 'The item you specified is not valid. Please specify one of the following variablesas a string:\n' \
                  + 'gene_type\n' \
                  + 'TF_list\n' \
                  + 'DG_universe\n' \
                  + 'DG_TF\n' \
                  + 'DEG_list\n' \
                  + 'DEG_filename\n' \
                  + 'DEG_to_pvalue\n' \
                  + 'DEG_to_updown\n' \
                  + 'DEG_full_graph\n' \
                  + 'tf_target_enrichment\n' \
                  + 'tf_enrichment\n' \
                  + 'z_scores\n\n'
            return None

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
        elif item == 'TF_list': self.TF_list = value
        elif item == 'DG_universe': self.DG_universe = value
        elif item == 'DG_TF': self.DG_TF = value
        elif item == 'DEG_list': self.DEG_list = value
        elif item == 'DEG_filename': self.DEG_filename = value
        elif item == 'DEG_to_pvalue': self.DEG_to_pvalue = value
        elif item == 'DEG_to_updown': self.DEG_to_updown = value
        elif item == 'DEG_full_graph': self.DEG_full_graph = value
        elif item == 'tf_target_enrichment': self.tf_target_enrichment = value
        elif item == 'tf_enrichment': self.tf_enrichment = value
        elif item == 'z_scores': self.z_scores = value
        else:
            print 'The item you specified is not valid. Please specify one of the following variables:\n' \
                  + 'gene_type\n' \
                  + 'TF_list\n' \
                  + 'DG_universe\n' \
                  + 'DG_TF\n' \
                  + 'DEG_list\n' \
                  + 'DEG_filename\n' \
                  + 'DEG_to_pvalue\n' \
                  + 'DEG_to_updown\n' \
                  + 'DEG_full_graph\n' \
                  + 'tf_target_enrichment\n' \
                  + 'tf_enrichment\n' \
                  + 'z_scores\n\n'




    # ----------------------------- TRANSCRIPTION FACTOR -------------------------- #

    def load_slowkow(self, filename_list=['../TF_databases/slowkow_databases/TRED_TF.txt',
                                '../TF_databases/slowkow_databases/ITFP_TF.txt',
                                '../TF_databases/slowkow_databases/ENCODE_TF.txt',
                                '../TF_databases/slowkow_databases/Neph2012_TF.txt',
                                '../TF_databases/slowkow_databases/TRRUST_TF.txt',
                                '../TF_databases/slowkow_databases/Marbach2016_TF.txt']):
        TF_list = create_graph.load_slowkow(filename_list)
        self.TF_list = TF_list

    def load_jaspar(self, filename='jaspar_genereg_matrix.txt'):
        TF_list = create_graph.load_jaspar(filename)
        self.TF_list = TF_list

    def easy_load_TF_list(self, slowkow_bool=True,
                          slowkow_files=['../TF_databases/slowkow_databases/TRED_TF.txt',
                                         '../TF_databases/slowkow_databases/ITFP_TF.txt',
                                         '../TF_databases/slowkow_databases/ENCODE_TF.txt',
                                         '../TF_databases/slowkow_databases/Neph2012_TF.txt',
                                         '../TF_databases/slowkow_databases/TRRUST_TF.txt',
                                         '../TF_databases/slowkow_databases/Marbach2016_TF.txt'],
                          jaspar_bool=True,
                          jaspar_file="../TF_databases/jaspar_genereg_matrix.txt"):

        TF_list = create_graph.easy_load_TF_list(slowkow_bool, slowkow_files, jaspar_bool, jaspar_file, self.gene_type)
        self.TF_list = TF_list


    # ------------------------- BACKGROUND NETWORK ------------------------------------ #


    def load_small_STRING_to_digraph(self, filename, TF_list=[]):
        DG_TF, DG_universe = create_graph.load_small_STRING_to_digraph(filename, TF_list)
        self.DG_TF = DG_TF
        self.DG_universe = DG_universe


    def load_STRING_to_digraph(self, filename, TF_list, confidence_filter=400):
        DG_TF, DG_universe = create_graph.load_STRING_to_digraph(filename, TF_list, confidence_filter, self.gene_type)
        self.DG_TF = DG_TF
        self.DG_universe = DG_universe




    # --------------------- DEG LOAD FUNCTIONS ---------------------------#


    def create_DEG_list(self, filename,

                        p_value_filter=0.05,
                        p_value_or_adj='adj',  # filtering by p-value ('p') or adjusted p-value ('adj')

                        fold_change_filter=None,  # specify a number to filter by absolute (log) fold change

                        gene_column_header=None,
                        p_value_column_header=None,
                        fold_change_column_header=None):

        # make sure user has run all prerequisites
        for item in [self.DG_TF]:
            if self.check_exists(item) == False:
                return -1

        # create the DEG list with specified cut-offs
        DEG_list, DG_TF = create_graph.create_DEG_list(filename, self.DG_TF, p_value_filter, p_value_or_adj,
                fold_change_filter, self.gene_type, gene_column_header, p_value_column_header, fold_change_column_header)

        self.DEG_list = DEG_list
        self.DEG_filename = filename
        self.DG_TF = DG_TF


        # create the full graph (call same function but just don't filter it)
        DEG_full_graph, DEG_to_pvalue, DEG_to_updown = create_graph.create_DEG_full_graph(filename,
                    p_value_or_adj=p_value_or_adj,
                    gene_symbol = self.gene_symbol,
                    gene_column_header=gene_column_header,
                    p_value_column_header=p_value_column_header,
                    fold_change_column_header=fold_change_column_header
                    )
        self.DEG_full_graph = DEG_full_graph
        self.DEG_to_pvalue = DEG_to_pvalue # keep track of all genes' (adj) p-values
        self.DEG_to_updown = DEG_to_updown # keep track of all genes' (log) fold change information



    # --------------------- P-VALUE FUNCTIONS --------------------------- #



    def tf_target_enrichment(self):

        # make sure user has run all prerequisites
        for item in [self.DG_TF, self.DG_universe, self.DEG_list]:
            if self.check_exists(item) == False:
                return

        tf_target_enrichment = stat_analysis.tf_target_enrichment(self.DG_TF, self.DG_universe, self.DEG_list)
        self.tf_target_enrichment = tf_target_enrichment


    def tf_enrichment(self):

        # make sure user has run all prerequisites
        for item in [self.TF_list, self.DG_full_graph, self.DEG_list]:
            if self.check_exists(item) == False:
                return

        tf_enrichment = stat_analysis.tf_enrichment(self.TF_list, self.DEG_full_graph, self.DEG_list)
        self.tf_enrichment = tf_enrichment



    # --------------------- Z-SCORE FUNCTIONS --------------------------- #

    def tf_zscore(self, bias_filter=0.25):

        # make sure user has run all prerequisites
        for item in [self.DG_TF, self.DEG_list]:
            if self.check_exists(item) == False:
                return

        z_scores = stat_analysis.tf_zscore(self.DG_TF, self.DEG_list, bias_filter)
        self.z_scores = z_scores


    # --------------------- DISPLAY FUNCTIONS --------------------------- #

    # TODO: fix this header and function to print adjp and fold change
    def top_values(self, act=True, abs_value=False, top=10):

        # make sure user has run all prerequisites
        for item in [self.z_scores, self.DEG_to_pvalue, self.DEG_to_updown]:
            if self.check_exists(item) == False:
                return -1

        return stat_analysis.top_values(self.z_scores, self.DEG_to_pvalue, self.DEG_to_updown, act, abs_value, top)

    # TODO: fix this make histogram
    def rank_and_score_df(self, series, genes_to_rank, value_name='z-score', abs_value=True, act=False, remove_dups=False):
        return stat_analysis.rank_and_score_df(series, genes_to_rank, value_name, abs_value, act, remove_dups)

    def vis_tf_network(self, tf,
                       directed_edges=False,
                       node_spacing=2200,
                       color_non_DEGs=False,
                       color_map=plt.cm.bwr,
                       graph_id=0
                       ):

        # make sure user has run all prerequisites
        for item in [self.DG_TF, self.DEG_filename, self.DEG_list]:
            if self.check_exists(item) == False:
                return -1

         return stat_analysis.vis_tf_network(self.DG_TF, tf, self.DEG_filename, self.DEG_list, directed_edges, node_spacing, color_non_DEGs, color_map, graph_id)

    # ------------------- HELPER FUNCTIONS ------------------------------ #


    def check_exists(self, item):

        item_to_message = {}

        item_to_message[self.TF_list] = 'No transcription factor list currently on file. Please run one of the following methods:\n'\
                                        + 'Upstream.load_slowkow\n'\
                                        + 'Upstream.load_jaspar\n'\
                                        + 'Upstream.easy_load_TF_list\n'\
                                        + 'Or assign your own using self.TF_list\n'
        item_to_message[self].DG_universe = 'No background network currently on file. Please run one of the following methods:\n'\
                                        + 'Upstream.load_small_STRING_to_digraph\n'\
                                        + 'Upstream.load_STRING_to_digraph\n'\
                                        + 'Or assign your own using self.DG_universe\n'
        item_to_message[self].DG_TF = 'No transcription factor background network currently on file. Please run one of the following methods:\n'\
                                        + 'Upstream.load_small_STRING_to_digraph\n'\
                                        + 'Upstream.load_STRING_to_digraph\n'\
                                        + 'Or assign your own using self.DG_TF\n'
        item_to_message[self].DEG_list = 'No differentially expressed gene list currently on file. Please run the following method:\n'\
                                        + 'Upstream.create_DEG_list\n'\
                                        + 'Or assign your own using self.DEG_list\n'
        item_to_message[self].DEG_filename = 'No differentially expressed gene filename currently on file. Please run the following method:\n' \
                                        + 'Upstream.create_DEG_list\n'\
                                        + 'Or assign your own using self.DEG_filename\n'
        item_to_message[self].DEG_to_pvalue = 'No DEG to (adj) p-value mapping currently on file. Please run the following method:\n'\
                                        + 'Upstream.create_DEG_list\n'\
                                        + 'Or assign your own using self.DEG_to_pvalue\n'
        item_to_message[self].DEG_to_updown = 'No DEG to up/down regulation mapping currently on file. Please run the following method:\n'\
                                        + 'Upstream.create_DEG_list\n'\
                                        + 'Or assign your own using self.DEG_to_updown\n'
        item_to_message[self].DEG_full_graph = 'No full expression gene graph currently on file. Please run the following method:\n'\
                                        + 'Upstream.create_DEG_list\n'\
                                        + 'Or assign your own using self.DEG_full_graph\n'
        item_to_message[self].tf_target_enrichment = 'No TF-target enrichment dictionary currently on file. Please run the following method:\n'\
                                        + 'Upstream.tf_target_enrichment\n'\
                                        + 'Or assign your own using self.tf_target_enrichment\n'
        item_to_message[self].tf_enrichment = 'No TF enrichment value currently on file. Please run the following method:\n'\
                                        + 'Upstream.tf_enrichment\n'\
                                        + 'Or assign your own using self.tf_enrichment\n'
        item_to_message[self].z_scores = 'No z-score dictionary currently on file. Please run the following method:\n'\
                                        + 'Upstream.tf_zscore\n'\
                                        + 'Or assign your own using self.z_score\n'
        if item == None:
            print item_to_message[item]
            return False

        return True