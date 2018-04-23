"""
-------------------------------------------
Author: Mikayla Webster (13webstermj@gmail.com)
Date: 2/5/18
-------------------------------------------
"""

import create_graph
reload(create_graph)
import stat_analysis
reload(stat_analysis)
import matplotlib.pyplot as plt
import copy

class Upstream:

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

        # create_graph variables
        self.TF_list = None
        self.DG_universe = None
        self.DG_TF = None
        self.DEG_list = None
        self.DEG_filename = None
        self.DEG_to_pvalue = None # (adj) p-value of all genes
        self.DEG_to_updown = None # (log) fold change of all genes
        self.DEG_full_graph = None

        # stat_analysis variables
        self.tf_target_enrichment = None
        self.tf_enrichment = None
        self.z_scores = None

        # map string to actual instance variable
        self.string_to_item = {}
        self.string_to_item['gene_type'] = self.gene_type
        self.string_to_item['species'] = self.species
        self.string_to_item['TF_list'] = self.TF_list
        self.string_to_item['DG_universe'] = self.DG_universe
        self.string_to_item['DG_TF'] = self.DG_TF
        self.string_to_item['DEG_list'] = self.DEG_list
        self.string_to_item['DEG_filename'] = self.DEG_filename
        self.string_to_item['DEG_to_pvalue'] = self.DEG_to_pvalue
        self.string_to_item['DEG_to_updown'] = self.DEG_to_updown
        self.string_to_item['DEG_full_graph'] = self.DEG_full_graph
        self.string_to_item['tf_target_enrichment'] = self.tf_target_enrichment
        self.string_to_item['tf_enrichment'] = self.tf_enrichment
        self.string_to_item['z_scores'] = self.z_scores

        # instanciating error message dict
        self.item_to_message = {}
        self.item_to_message['gene_type'] = 'No gene type currently on file. Please create a new instance of Upstream using '\
                             + 'Upstream(gene_type, species)'
        self.item_to_message['species'] = 'No species declaration currently on file. Please create a new instance of Upstream using ' \
                             + 'Upstream(gene_type, species)'
        self.item_to_message['TF_list'] = 'No transcription factor list currently on file. Please run one of the following methods:\n' \
                             + ' - Upstream.load_slowkow\n' \
                             + ' - Upstream.load_jaspar\n' \
                             + ' - Upstream.easy_load_TF_list\n' \
                             + 'Or assign your own using self.TF_list\n'
        self.item_to_message['DG_universe'] = 'No background network currently on file. Please run one of the following methods:\n' \
                             + ' - Upstream.load_small_STRING_to_digraph\n' \
                             + ' - Upstream.load_STRING_to_digraph\n' \
                             + 'Or assign your own using self.DG_universe\n'
        self.item_to_message['DG_TF'] = 'No transcription factor background network currently on file. Please run one of the following methods:\n' \
                             + ' - Upstream.load_small_STRING_to_digraph\n' \
                             + ' - Upstream.load_STRING_to_digraph\n' \
                             + 'Or assign your own using self.DG_TF\n'
        self.item_to_message['DEG_list'] = 'No differentially expressed gene list currently on file. Please run the following method:\n' \
                             + ' - Upstream.create_DEG_list\n' \
                             + 'Or assign your own using self.DEG_list\n'
        self.item_to_message['DEG_filename'] = 'No differentially expressed gene filename currently on file. Please run the following method:\n' \
                             + ' - Upstream.create_DEG_list\n' \
                             + 'Or assign your own using self.DEG_filename\n'
        self.item_to_message['DEG_to_pvalue'] = 'No DEG to (adj) p-value mapping currently on file. Please run the following method:\n' \
                             + ' - Upstream.create_DEG_list\n' \
                             + 'Or assign your own using self.DEG_to_pvalue\n'
        self.item_to_message['DEG_to_updown'] = 'No DEG to up/down regulation mapping currently on file. Please run the following method:\n' \
                             + ' - Upstream.create_DEG_list\n' \
                             + 'Or assign your own using self.DEG_to_updown\n'
        self.item_to_message['DEG_full_graph'] = 'No full expression gene graph currently on file. Please run the following method:\n' \
                             + ' - Upstream.create_DEG_list\n' \
                             + 'Or assign your own using self.DEG_full_graph\n'
        self.item_to_message['tf_target_enrichment'] = 'No TF-target enrichment dictionary currently on file. Please run the following method:\n' \
                             + ' - Upstream.tf_target_enrichment\n' \
                             + 'Or assign your own using self.tf_target_enrichment\n'
        self.item_to_message['tf_enrichment'] = 'No TF enrichment value currently on file. Please run the following method:\n' \
                             + ' - Upstream.tf_enrichment\n' \
                             + 'Or assign your own using self.tf_enrichment\n'
        self.item_to_message['z_scores'] = 'No z-score dictionary currently on file. Please run the following method:\n' \
                             + ' - Upstream.tf_zscore\n' \
                             + 'Or assign your own using self.z_score\n'

    def copyUpstream(self):

        to_return = Upstream(self.gene_type, self.species)

        # create_graph variables
        to_return.TF_list = self.TF_list
        to_return.DG_universe = copy.deepcopy(self.DG_universe)
        to_return.DG_TF = copy.deepcopy(self.DG_TF)
        to_return.DEG_list = self.DEG_list
        to_return.DEG_filename = self.DEG_filename
        to_return.DEG_to_pvalue = self.DEG_to_pvalue
        to_return.DEG_to_updown = self.DEG_to_updown
        to_return.DEG_full_graph = copy.deepcopy(self.DEG_full_graph)

        # stat_analysis variables
        to_return.tf_target_enrichment = copy.deepcopy(self.tf_target_enrichment)
        to_return.tf_enrichment = copy.deepcopy(self.tf_enrichment)
        to_return.z_scores = copy.deepcopy(self.z_scores)

        # misc
        to_return.string_to_item = self.string_to_item
        to_return.item_to_message = self.item_to_message

        return to_return


    def __repr__(self):
        return self.create_to_string()

    def __str__(self):
        return self.create_to_string()

    def create_to_string(self):
        for item in ['gene_type', 'species', 'TF_list', 'DG_universe', 'DG_TF', 'DEG_list', 'DEG_filename',
                     'DEG_to_pvalue', 'DEG_to_updown', 'DEG_full_graph', 'tf_target_enrichment',
                     'tf_enrichment', 'z_scores']:
            print item + ': '
            exists = self.check_exists(item)
            if exists == True:
                print str(type(self.string_to_item[item])) + '\n'
        return '\n'


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
            print 'The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
                  + 'gene_type\n' \
                  + 'species\n' \
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

    def easy_load_TF_list(self, csv_filename, jaspar = True, TRED = True, ITFP = True, ENCODE = True, 
                          Neph2012 = True, TRRUST = True, Marbach2016 = True):

        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        TF_list = create_graph.easy_load_TF_list(csv_filename, jaspar, TRED, ITFP, ENCODE, Neph2012, TRRUST, Marbach2016, self.species, self.gene_type)
        self.TF_list = TF_list


    # ------------------------- BACKGROUND NETWORK ------------------------------------ #


    def load_small_STRING_to_digraph(self, filename): # TODO: check if we can entrez/species

        # make sure user has run all prerequisites
        if self.check_exists('TF_list') == False:
            return

        DG_TF, DG_universe = create_graph.load_small_STRING_to_digraph(filename, self.TF_list)
        self.DG_TF = DG_TF
        self.DG_universe = DG_universe


    def load_STRING_to_digraph(self, filename, confidence_filter=400):

        # make sure user has run all prerequisites
        for item in ['TF_list', 'gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        DG_TF, DG_universe = create_graph.load_STRING_to_digraph(filename, self.TF_list, confidence_filter, self.gene_type, self.species)
        self.DG_TF = DG_TF
        self.DG_universe = DG_universe




    # --------------------- DEG LOAD FUNCTIONS ---------------------------#


    def create_DEG_list(self, filename,

                        p_value_filter=0.05,
                        p_value_or_adj='adj',  # filtering by p-value ('p') or adjusted p-value ('adj')
                        fold_change_filter=None,  # specify a number to filter by absolute (log) fold change

                        gene_column_header=None,
                        p_value_column_header=None,
                        fold_change_column_header=None,
						sep = '\t'):

        # make sure user has run all prerequisites
        for item in ['DG_TF', 'gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        # create the DEG list with specified cut-offs
        DEG_list, DG_TF = create_graph.create_DEG_list(filename, self.DG_TF, p_value_filter, p_value_or_adj,
                                                       fold_change_filter, self.gene_type, gene_column_header, p_value_column_header, fold_change_column_header, sep)

        self.DEG_list = DEG_list
        self.DEG_filename = filename
        self.DG_TF = DG_TF


        # create the full graph (call same function but just don't filter it)
        DEG_full_graph, DEG_to_pvalue, DEG_to_updown = create_graph.create_DEG_full_graph(filename,
                                                                                          p_value_or_adj=p_value_or_adj,
                                                                                          gene_type = self.gene_type,
                                                                                          gene_column_header=gene_column_header,
                                                                                          p_value_column_header=p_value_column_header,
                                                                                          fold_change_column_header=fold_change_column_header,
                                                                                          sep = sep
                                                                                          )
        self.DEG_full_graph = DEG_full_graph
        self.DEG_to_pvalue = DEG_to_pvalue # keep track of all genes' (adj) p-values
        self.DEG_to_updown = DEG_to_updown # keep track of all genes' (log) fold change information



    # --------------------- P-VALUE FUNCTIONS --------------------------- #



    def tf_target_enrichment_calc(self):

        # make sure user has run all prerequisites
        for item in ['DG_TF', 'DG_universe', 'DEG_list']:
            if self.check_exists(item) == False:
                return

        tf_target_enrichment = stat_analysis.tf_target_enrichment(self.DG_TF, self.DG_universe, self.DEG_list)
        self.tf_target_enrichment = tf_target_enrichment


    def tf_enrichment_calc(self):

        # make sure user has run all prerequisites
        for item in ['TF_list', 'DEG_full_graph', 'DEG_list']:
            if self.check_exists(item) == False:
                return

        tf_enrichment = stat_analysis.tf_enrichment(self.TF_list, self.DEG_full_graph, self.DEG_list)
        self.tf_enrichment = tf_enrichment



    # --------------------- Z-SCORE FUNCTIONS --------------------------- #

    def tf_zscore(self, bias_filter=0.25):

        # make sure user has run all prerequisites
        for item in ['DG_TF', 'DEG_list']:
            if self.check_exists(item) == False:
                return

        z_scores = stat_analysis.tf_zscore(self.DG_TF, self.DEG_list, bias_filter)
        self.z_scores = z_scores


    # --------------------- DISPLAY FUNCTIONS --------------------------- #

    def top_values(self, act=True, abs_value=False, top=10):

        # make sure user has run all prerequisites
        for item in ['z_scores', 'DEG_to_pvalue', 'DEG_to_updown']:
            if self.check_exists(item) == False:
                return -1

        return stat_analysis.top_values(self.z_scores, self.DEG_to_pvalue, self.DEG_to_updown, act, abs_value, top)


    def compare_genes(self, genes_to_rank, fig_size=(12, 7), font_size=12, anno_vert_dist=0.025):

        # make sure user has run all prerequisites
        if self.check_exists('z_scores') == False:
            return -1

        stat_analysis.compare_genes(self.z_scores, genes_to_rank, fig_size=fig_size, font_size=font_size, anno_vert_dist=anno_vert_dist)


    def vis_tf_network(self, tf,
                       directed_edges=False,
                       node_spacing=2200,
                       color_non_DEGs=False,
                       color_map=plt.cm.bwr,
                       graph_id=0,
                       tf_size_amplifier = 8,
					   alpha = 1.0, 
				       color_vals_transform = None,
				       ceil_val=10,
                       color_max_frac = 1.0,
				       color_min_frac = 0.0,
				       vmin=None,
				       vmax=None,
					   tf_shape = 'star'
                       ):

        # make sure user has run all prerequisites
        for item in ['DG_TF', 'DEG_list', 'DEG_to_pvalue', 'DEG_to_updown']:
            if self.check_exists(item) == False:
                return -1

        return stat_analysis.vis_tf_network(self.DG_TF, tf, self.DEG_list, self.DEG_to_pvalue, self.DEG_to_updown, directed_edges, node_spacing, color_non_DEGs, color_map, graph_id, tf_size_amplifier, alpha = alpha, color_vals_transform = color_vals_transform,
                                            ceil_val = ceil_val, color_max_frac = color_max_frac, color_min_frac = color_min_frac, vmin = vmin, vmax = vmax, tf_shape = tf_shape)




    def to_csv(self, out_filename):

        # make sure user has run all prerequisites
        for item in ['z_scores', 'DEG_to_pvalue', 'DEG_to_updown', 'tf_target_enrichment', 'DG_TF']:
            if self.check_exists(item) == False:
                return -1

        stat_analysis.to_csv(out_filename, self.z_scores, self.DEG_to_pvalue, self.DEG_to_updown, self.tf_target_enrichment, self.DG_TF)


    # ------------------- HELPER FUNCTIONS ------------------------------ #


    # item must not be strign version
    def check_exists(self, item):

        # re-map it so it stays up to date
        self.string_to_item['gene_type'] = self.gene_type
        self.string_to_item['species'] = self.species
        self.string_to_item['TF_list'] = self.TF_list
        self.string_to_item['DG_universe'] = self.DG_universe
        self.string_to_item['DG_TF'] = self.DG_TF
        self.string_to_item['DEG_list'] = self.DEG_list
        self.string_to_item['DEG_filename'] = self.DEG_filename
        self.string_to_item['DEG_to_pvalue'] = self.DEG_to_pvalue
        self.string_to_item['DEG_to_updown'] = self.DEG_to_updown
        self.string_to_item['DEG_full_graph'] = self.DEG_full_graph
        self.string_to_item['tf_target_enrichment'] = self.tf_target_enrichment
        self.string_to_item['tf_enrichment'] = self.tf_enrichment
        self.string_to_item['z_scores'] = self.z_scores

        try:
            if (type(self.string_to_item[item]) == type(None)):
                print self.item_to_message[item]
                return False
        except:
            print 'The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
                  + '- gene_type\n' \
                  + '- species\n' \
                  + '- TF_list\n' \
                  + '- DG_universe\n' \
                  + '- DG_TF\n' \
                  + '- DEG_list\n' \
                  + '- DEG_filename\n' \
                  + '- DEG_to_pvalue\n' \
                  + '- DEG_to_updown\n' \
                  + '- DEG_full_graph\n' \
                  + '- tf_target_enrichment\n' \
                  + '- tf_enrichment\n' \
                  + '- z_scores\n\n'
            return False

        return True