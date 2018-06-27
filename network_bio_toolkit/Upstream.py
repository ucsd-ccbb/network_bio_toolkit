"""
-------------------------------------------
Author: Mikayla Webster (13webstermj@gmail.com)
Date: 2/5/18
-------------------------------------------
"""

# our modules
import create_graph
import stat_analysis

# common packages, most likely already installed
import matplotlib.pyplot as plt
import copy

# reloading, for testing
reload(create_graph)
reload(stat_analysis)

class Upstream:

    def __init__(self, gene_type = 'symbol', species = 'human'):
    
        """
            Initiates Upstream class instance, initiating all instance variables to None.
        """

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
    
        """
            Creates a deep copy of this analysis by initiating a new Upstream instance
            and (deep) copying all instance variables to this new instance.
        """

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
    
        """
            Creates a string representation of this object by printing whether each of its
            instance variables has been instanciated or not. If an instance variable has been 
            instanciated, it will print the type of that variable.
        """
    
        for item in ['gene_type', 'species', 'TF_list', 'DG_universe', 'DG_TF', 'DEG_list', 'DEG_filename',
                     'DEG_to_pvalue', 'DEG_to_updown', 'DEG_full_graph', 'tf_target_enrichment',
                     'tf_enrichment', 'z_scores']:
            print(item + ': ')
            exists = self.check_exists(item)
            if exists == True:
                print(str(type(self.string_to_item[item])) + '\n')
        return '\n'


    # ----------------------------- GETTERS AND SETTERS -------------------------- #
    # So our users don't accidentally modify an instance variable without realizing

    # to ask to see a variable
    def get(self, item):
    
        """
            Use this function to access instance variable values, but not the actual instance variable
            references. (aka, if you modify when this function returns, it will not modify the stored
            version of that variable.)
            
            Args:
                item: String, name of the variable you'd like to access. (Look at ura_instance.item_to_message 
                      if you are unsure.
        """

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
    
        """
            Use this function to update the value of an instance variable.
            
            Args:
                item: String, name of the variable you'd like to modify. (Look at ura_instance.item_to_message 
                      if you are unsure.
        """

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
            print('The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
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
                  + 'z_scores\n\n')




    # ----------------------------- TRANSCRIPTION FACTOR -------------------------- #

    def easy_load_TF_list(self, csv_filename, jaspar = True, TRED = True, ITFP = True, ENCODE = True, 
                          Neph2012 = True, TRRUST = True, Marbach2016 = True):
                          
        """
            Loads a list of transcription factors (TF) containing genes from the user specified databases.

            Args:
                csv_filename: String, filepath of where to find TF data ('../../TF_databases/TF_database_URA.csv' if you are using our 
                    github directory structure)
                jaspar: Boolean, whether or not to include jaspar database in TF list
                TRED: Boolean, whether or not to include TRED database in TF list
                ITFP: Boolean, whether or not to include ITFP database in TF list
                ENCODE: Boolean, whether or not to include ENCODE database in TF list
                Neph2012: Boolean, whether or not to include Neph2012 database in TF list
                TRRUST: Boolean, whether or not to include TRRUST database in TF list
                Marbach2016: Boolean, whether or not to include Marbach2016 database in TF list
        """

        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        TF_list = create_graph.easy_load_TF_list(csv_filename, jaspar, TRED, ITFP, ENCODE, Neph2012, TRRUST, Marbach2016, self.species, self.gene_type)
        self.TF_list = TF_list


    # ------------------------- BACKGROUND NETWORK ------------------------------------ #


    def load_small_STRING_to_digraph(self, filename): # TODO: check if we can entrez/species
    
        """
            ** This is a depricated function that may not be compatible with the rest of our functions
        
            This function loads a small subset of the STRING database into a networkx digraph that can be used as input
            into the rest of our functions. This function filters the input database so that only the genes
            indicated by TF_list, and all of their out-going neighbors, will remain in the graph. Namely, the only
            source nodes left in the graph will be the genes in TF_list.

            ** Note that there is a loss of 110 our of about 40000 edges due to removal of multiedges **

            Args:
                filename: String, the path to the string file
        """

        # make sure user has run all prerequisites
        if self.check_exists('TF_list') == False:
            return

        DG_TF, DG_universe = create_graph.load_small_STRING_to_digraph(filename, self.TF_list)
        self.DG_TF = DG_TF
        self.DG_universe = DG_universe


    def load_STRING_to_digraph(self, filename, confidence_filter = 400):
    
        """
            This function loads a subset of the STRING database from either '9606.protein.actions.v10.5.txt' or
            '10090.protein.actions.v10.5.txt' into a networkx digraph that can be used as the background network for 
            the rest of our functions. IT ONLY KEEPS ACTIVATING AND INHIBITING EDGES FROM THE STRING NETWORK. This 
            function filters the input database so that only the genes indicated by TF_list, and all of their 
            out-going neighbors, will remain in the graph. Namely, the only source nodes left in the graph will be the 
            genes in TF_list.

            Args:
                filename: String, the filepath to the STRING database file.
                confidence_filter: A number between 0 and 1000, all interactions with confidece less than this number will be filtered out
        """

        # make sure user has run all prerequisites
        for item in ['TF_list', 'gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        DG_TF, DG_universe = create_graph.load_STRING_to_digraph(filename, self.TF_list, confidence_filter, self.gene_type, self.species)
        self.DG_TF = DG_TF
        self.DG_universe = DG_universe
        
    
    def load_ndex_from_server(self, UUID, relabel_node_field = None):
    
        """
            This function loads an Ndex network as a networkx graph. This function can also filter the input database 
            so that only the genes indicated by filter_list, and all of their out-going neighbors, will remain in the 
            graph. This function is comparible to load_STRING_to_digraph().

            Args:
                UUID: String, the unique ID associated with the Ndex network you wish to load.
                      Example: 'e11c6684-5ac2-11e8-a4bf-0ac135e8bacf' <-- STRING_human_protein_actions_confidence_700_edges network
                               '76160faa-5d0f-11e8-a4bf-0ac135e8bacf' <-- STRING_mouse_protein_links_confidence_700_edges network
                               
                relabel_node_field: String, node attribute to relable the nodes in this graph with. Most Ndex graphs have a 'name' field
                      that contains the gene name you are looking for for each node. (Sometimes the nodes themselves are named weirdly).
                      In that case, set relabel_node_field = 'name'.
        """
    
        # make sure user has run all prerequisites
        for item in ['gene_type', 'species']:
            if self.check_exists(item) == False:
                return
    
        DG_universe = create_graph.load_ndex_from_server(UUID, relabel_node_field, None)
        self.DG_universe = DG_universe




    # --------------------- DEG LOAD FUNCTIONS ---------------------------#


    def create_DEG_list(self, 
                        filename,
                        p_value_filter = 0.05,
                        p_value_or_adj = 'adj',  # filtering by p-value ('p') or adjusted p-value ('adj')
                        fold_change_filter = None,  # specify a number to filter by absolute (log) fold change
                        gene_column_header = None,
                        p_value_column_header = None,
                        fold_change_column_header = None,
						sep = '\t',
                        return_full_values = False):
                        
        """
            This function takes a differential expression data file and loads it into multiple dictionaries 
            that can be used by our other functions. Differential expression file must contain a gene name column, 
            a (log) fold change column, and an (adjusted) p-value column, each specified with a common column header.
            If your column header for any of these columns is not intuitive, specify the name of that column header.
            Each line in your input file should represent a different gene.

            Args:
                filename: String, filepath to the the differential expression input file
                p_value_filter: Float, typically a number between 0 and 1, the number to filter the adjusted p-vlaue by. Will remove all above this threshold
                p_value_or_adj: String, either 'p' or 'adj', that specifies whether we should use p-value or adjusted p-value information
                fold_change_filter: Float/Int, will filter out genes with fold change absolute value greater than this number
                gene_column_header: String, if your gene name column header is not intuitive, specify it here
                p_value_column_header: String, if your p-value column header is not intuitive, specify it here
                fold_change_column_header: String, if your fold change column header is not intuitive, specify it here
                sep: String, separating agent in your input file. Common exmaples include '\t' and ','
                return_full_values: Boolean, specifies whether to return only Differentially Expressed Genes' p-value and log fold
                    information, or whether to return entire file's information.
        """

        # make sure user has run all prerequisites
        for item in ['DG_TF', 'gene_type', 'species']:
            if self.check_exists(item) == False:
                return

        # create the DEG list with specified cut-offs
        DEG_list, DG_TF = create_graph.create_DEG_list(filename, self.DG_TF, p_value_filter, p_value_or_adj,
                                                       fold_change_filter, self.gene_type, gene_column_header, 
                                                       p_value_column_header, fold_change_column_header, sep,
                                                       return_full_values)

        self.DEG_list = DEG_list
        self.DEG_filename = filename
        self.DG_TF = DG_TF


        # create the full graph (call same function but just don't filter it)
        DEG_full_graph, DEG_to_pvalue, DEG_to_updown = create_graph.create_DEG_full_graph(filename,
                                                                                          p_value_or_adj = p_value_or_adj,
                                                                                          gene_type = self.gene_type,
                                                                                          gene_column_header = gene_column_header,
                                                                                          p_value_column_header = p_value_column_header,
                                                                                          fold_change_column_header = fold_change_column_header,
                                                                                          sep = sep
                                                                                          )
        self.DEG_full_graph = DEG_full_graph
        self.DEG_to_pvalue = DEG_to_pvalue # keep track of all genes' (adj) p-values
        self.DEG_to_updown = DEG_to_updown # keep track of all genes' (log) fold change information

        
        
    # ----------------------- LOCALIZATION -------------------------------#
    
    def localization(self, num_reps = 10, sample_frac = 0.8, method = 'numedges', plot = True, print_counter = False):
    
        """
            Function to calculate localization of your DEGs on your background network.
            Option to compute number of edges (method = 'numedges') or largest connected component (method = 'LLC') 
            localization analysis. Calculates by sampling sub-sections of the focal genes/random set. Percentage to sample
            is set by sample_frac. Option to plot the distributions of random and focal gene localizaiton.
            
            Args:
                num_reps: Int, number of times to randomly sample
                sample_frac: Float, percent of sampled genes
                method: String, to decide which type of localization analysis to run. Options: 'numedges', 'LLC', or 'both'.
                plot: Bool, whether to plot the distributions in the output jupyter notebook cell
                print_counter: Bool, whether to print a counter that tells you which iteration you are on (every 25 iterations).
                               Useful when the num_reps is very high.
                
            Returns: 
                numedges_list: List, the number of edges calculated for each rep, sampling over focal genes. 
                    Empty if method = 'LLC'. 
                numedges_rand: List, the number of edges calculated for each rep, sampling over random genes of 
                    similar degree in the background network. Empty if method = 'LLC'.
                LCC_list: List, the size of the largest connected component, calculated for each rep, sampling over focal genes. 
                    Empty if method = 'numedges'. 
                LCC_rand: List, the size of the largest connected component, calculated for each rep, sampling over random genes of 
                    similar degree in the background network. Empty if method = 'numedges'. 
        """
    
        # make sure user has run all prerequisites
        for item in ['DG_universe', 'DEG_list']:
            if self.check_exists(item) == False:
                return
    
        return stat_analysis.localization(self.DG_universe, self.DEG_list, num_reps, sample_frac, method, plot, print_counter)
        
    def localization_full(self, num_reps = 200, 
                          method = 'LCC', 
                          print_counter = False, 
                          label = 'focal genes',
                          line_height = 0.1,
                          legend_loc = 'upper left'):
    
        """
            Function to calculate localization of an input set of genes (focal_genes) on a background network (Gint).
            Option to compute number of edges (method = 'numedges') or largest connected component (method = 'LLC') 
            localization analysis. DOes no sub-sampling. Plots the distribution of random gene localizaiton, and 
            marks the focal set localization on distribution. Includes p-value of focal set localization.
            
            Args:
                num_reps: Int, number of times to randomly sample
                method: String, to decide which type of localization analysis to run. Options: 'numedges', 'LLC', or 'both'.
                print_counter: Bool, whether to print a counter that tells you which iteration you are on (every 25 iterations).
                               Useful when the num_reps is very high.
                label: String, label for focal genes in graph legend
                line_height: Float, the height of the red line that marks the focal gene localization
                legend_loc: String, relative position of legend in graph. Something similar to 'upper left'.
                
            Returns: 
                numedges_list: List, the number of edges calculated for each rep, over focal genes. 
                    Empty if method = 'LLC'. 
                numedges_rand: List, the number of edges calculated for each rep, over random genes of 
                    similar degree in the background network. Empty if method = 'LLC'.
                LCC_list: List, the size of the largest connected component, calculated for each rep, over focal genes. 
                    Empty if method = 'numedges'. 
                LCC_rand: List, the size of the largest connected component, calculated for each rep, over random genes of 
                    similar degree in the background network. Empty if method = 'numedges'. 
        """
    
        # make sure user has run all prerequisites
        for item in ['DG_universe', 'DEG_list']:
            if self.check_exists(item) == False:
                return
    
        return stat_analysis.localization_full(self.DG_universe, self.DEG_list, num_reps, method, print_counter, label, line_height, legend_loc)
        


    # --------------------- P-VALUE FUNCTIONS --------------------------- #



    def tf_target_enrichment_calc(self):
    
        """
            Our p-value function calculates the log of the p-value for every TF in the graph using [scipy.stats.hypergeom.logsf]
            (https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.hypergeom.html). These values help us
            determine which TF's are actually associated with our DEG's. If a TF is given a high value (because we are
            working with logs, not straight p-values), then it is likely that there is correlation between that TF and its
            DEG targets. Therefore, it is likely that TF is responsible for some of our observed gene expression.
            Note that if a TF is given a value of zero, that means none of the TF's targets were DEG's.
        """

        # make sure user has run all prerequisites
        for item in ['DG_TF', 'DG_universe', 'DEG_list']:
            if self.check_exists(item) == False:
                return

        tf_target_enrichment = stat_analysis.tf_target_enrichment(self.DG_TF, self.DG_universe, self.DEG_list)
        self.tf_target_enrichment = tf_target_enrichment


    def tf_enrichment_calc(self):
    
        """
            Our p-value function calculates the log of the p-value [scipy.stats.hypergeom.logsf] for all TF's in this analysis.
            (https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.hypergeom.html). This function generates
            a single value that represents how enriched the set of TF's is. If the p-value is high, then the set of TF's is enriched.
            (We are working with logs, not straight p-values)
        """

        # make sure user has run all prerequisites
        for item in ['TF_list', 'DEG_full_graph', 'DEG_list']:
            if self.check_exists(item) == False:
                return

        tf_enrichment = stat_analysis.tf_enrichment(self.TF_list, self.DEG_full_graph, self.DEG_list)
        self.tf_enrichment = tf_enrichment



    # --------------------- Z-SCORE FUNCTIONS --------------------------- #

    def tf_zscore(self, bias_filter = 0.25):
    
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
                bias_filter: number between 0 and 1, threshold to calculate z-score using biased formula
        """

        # make sure user has run all prerequisites
        for item in ['DG_TF', 'DEG_list']:
            if self.check_exists(item) == False:
                return

        z_scores = stat_analysis.tf_zscore(self.DG_TF, self.DEG_list, bias_filter)
        self.z_scores = z_scores


    # --------------------- DISPLAY FUNCTIONS --------------------------- #

    def top_values(self, act = True, abs_value = False, top = 10):
    
        """
            This function returns a sorted Pandas Series of the top (number indicated by top) genes based off z-score
            (given by z_score_series).

            Args:
                act: Boolean, True to sort by most positive z-score, False to sort by most negative. Ignored if abs_value
                    is True
                abs_value: Boolean, True to sort by absolute z-score, False otherwise
                top: the number of genes you wish to have returned

            Returns: A sorted Pandas Dataframe of the top genes, mapping each gene to its z-score
        """

        # make sure user has run all prerequisites
        for item in ['z_scores', 'DEG_to_pvalue', 'DEG_to_updown']:
            if self.check_exists(item) == False:
                return -1

        return stat_analysis.top_values(self.z_scores, self.DEG_to_pvalue, self.DEG_to_updown, act, abs_value, top)


    def compare_genes(self, genes_to_rank, fig_size = (12,7), font_size = 12, anno_vert_dist = 0.025):
    
        """
            Plots a histogram of the distribution of z-scores, and marks the genes indicated in genes_to_rank
            so that their values can be visually compared.

            Args:
                genes_to_rank: List, genes you wish to highlight in the distribution
                fig_size: (Int,Int), the dimensions of the histogram.
                font_size: Int, the font size of the genes_to_rank labels
                anno_vert_dist: Float, the distance separating the genes_to_rank labels from the x-axis

            Returns: A seaborn based histogram.
        """

        # make sure user has run all prerequisites
        if self.check_exists('z_scores') == False:
            return -1

        stat_analysis.compare_genes(self.z_scores, genes_to_rank, fig_size = fig_size, font_size = font_size, anno_vert_dist = anno_vert_dist)


    def vis_tf_network(self, tf,
                       directed_edges = False,
                       node_spacing = 2200,
                       color_non_DEGs = False,
                       color_map = plt.cm.bwr,
                       graph_id = 0,
                       tf_size_amplifier = 8,
					   alpha = 1.0, 
				       color_vals_transform = None,
				       ceil_val = 10,
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
                tf: String, all caps gene symbol of the regulator whose sub-network you wish to visualizes
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

        # make sure user has run all prerequisites
        for item in ['DG_TF', 'DEG_list', 'DEG_to_pvalue', 'DEG_to_updown']:
            if self.check_exists(item) == False:
                return -1

        return stat_analysis.vis_tf_network(self.DG_TF, tf, self.DEG_list, self.DEG_to_pvalue, self.DEG_to_updown, 
                                            directed_edges, node_spacing, color_non_DEGs, color_map, graph_id, 
                                            tf_size_amplifier, alpha = alpha, color_vals_transform = color_vals_transform,
                                            ceil_val = ceil_val, color_max_frac = color_max_frac, 
                                            color_min_frac = color_min_frac, vmin = vmin, vmax = vmax, tf_shape = tf_shape)




    def to_csv(self, out_filename):
    
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

        # make sure user has run all prerequisites
        for item in ['z_scores', 'DEG_to_pvalue', 'DEG_to_updown', 'tf_target_enrichment', 'DG_TF']:
            if self.check_exists(item) == False:
                return -1

        stat_analysis.to_csv(out_filename, self.z_scores, self.DEG_to_pvalue, self.DEG_to_updown, self.tf_target_enrichment, self.DG_TF)


    # ------------------- HELPER FUNCTIONS ------------------------------ #


    # item must not be strign version
    def check_exists(self, item):
    
        """
            Helper functions used to help guide the user as to which order they should call the analysis functions.
            Checks whether the named instance variable has been set yet.

            Args:
                item: String, name of the instance variable whose existance we are checking.
        """

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
                print(self.item_to_message[item])
                return False
        except:
            print('The item you specified (' + str(item) + ') is not valid. Please specify one of the following variables:\n' \
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
                  + '- z_scores\n\n')
            return False

        return True