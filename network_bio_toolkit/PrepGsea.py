import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
from ipywidgets import widgets
from IPython.display import display

class PrepGsea:
    
    def __init__(self, gmt_file, expression_file, meta_file, output_dir, display_df_meta = True):
        
        self.gmt_file = gmt_file 
        self.expression_file = expression_file
        self.meta_file = meta_file
        self.output_dir = output_dir
        
        self.df_meta = pd.read_csv(meta_file)
        if display_df_meta == True: # recommended so that you can see your class and control choices
            print "Meta DataFrame: "
            display(self.df_meta)
    
    def choose_class_col(self, sample_col_name = 'sample'):
        
        self.col_headers = list(self.df_meta)
        self.sample_col_name = sample_col_name
        
        self.col_headers = [header for header in self.col_headers if sample_col_name not in header.lower()]

        print 'Please choose a class from the following list.\n\n' \
                     + 'Do not reload cell after marking your choice.\n' \
                     + 'Mark your choice before proceding.'
        
        self.class_col_button = widgets.RadioButtons(
            options = self.col_headers,
            description='Class options:',
            style = {'description_width': 'initial'},
            disabled=False
        )
        
        display(self.class_col_button)
        
    def choose_two_classes(self):
        
        print 'Please specify the two classes you wish to compare for class \"' + str(self.class_col_button.value) + '\".\n\n' \
                     + 'Do not reload cell after typing your choices.\n' \
                     + 'Type your choices before proceding.\n' \
                     + 'Reload this cell if the correct class isn\'t showing'

        class_name = self.class_col_button.value
        self.class_A_textbox = widgets.Text(description = class_name + ' class 1: ', style = {'description_width': 'initial'})
        self.class_B_textbox = widgets.Text(description = class_name + ' class 2: ', style = {'description_width': 'initial'})
        display(self.class_A_textbox)
        display(self.class_B_textbox)
        
    def choose_controls(self):
        
        class_A = self.class_A_textbox.value
        class_B = self.class_B_textbox.value
        class_name = self.class_col_button.value

        if class_name in self.col_headers:
            self.col_headers.remove(class_name)

        print 'Please specify your control values for each of the following.\n\n' \
             + 'Do not reload cell after typing your choices.\n' \
             + 'Type your choices before proceding.'
        
        self.controls_widget_array = []
        for header in self.col_headers:
            if self.sample_col_name not in header.lower():
                text = widgets.Text(description = header, style = {'description_width': 'initial'})
                display(text)
                self.controls_widget_array.append(text)
                
    def filter_metafile(self):
        
        # remember user preferences
        class_A = self.class_A_textbox.value
        class_B = self.class_B_textbox.value
        class_name = self.class_col_button.value
                
        # capture control values
        self.controls = {}
        for col, text in zip(self.col_headers, self.controls_widget_array):
            self.controls[col] = text.value
            
        print 'Filtering metafile by the following classes and controls:\n\n' \
        + 'class name: ' + str(class_name) + '\n' \
        + str(class_name) + ' value 1: ' + str(class_A) + '\n' \
        + str(class_name) + ' value 2: ' + str(class_B) + '\n\n' \
        + 'controls: '
        
        for k,v in self.controls.items():
            print str(k) + ': ' + str(v)
        
        print '\nPlease confirm that all the information above is correct.\n'
        
        # filter by control values
        for key, value in self.controls.iteritems(): 
            self.df_meta = self.df_meta[self.df_meta[key]==value]
        
        # filter by class name
        self.samp_to_class = self.df_meta[['Sample_name', class_name]]
        self.samp_to_class = self.samp_to_class[(self.samp_to_class[class_name] == class_A) | (self.samp_to_class[class_name] == class_B)]
        
        display(self.df_meta)
        
    def filter_expression_file(self, sample_col_name = 'Sample_name'):
        
            # remember user preferences
            class_A = self.class_A_textbox.value
            class_B = self.class_B_textbox.value
            class_name = self.class_col_button.value

            # load and display expression file
            self.df_expression = pd.read_table(self.expression_file, index_col='Unnamed: 0')
            focal_samples = list(self.df_expression)  # header

            # filter netafile by samples
            self.df_meta = df_meta[df_meta[sample_col_name].isin(focal_samples)]

            # filter metafile by controls
            for key, value in controls.iteritems(): 
                self.df_meta = self.df_meta[self.df_meta[key]==value]

            print "Expression file before filtering: " + str(self.df_expression.shape)
            display(self.df_expression.head())
            print "\nFiltered meta file: " + str(self.df_meta.shape)
            display(self.df_meta)
            
            # extract only the COLUMNS with sample_name and class name
            self.samp_to_class = self.df_meta[['Sample_name', class_name]]

            # only keep the rows with class_A and class_B
            self.samp_to_class = self.samp_to_class[(self.samp_to_class[class_name] == class_A) | (self.samp_to_class[class_name] == class_B)]
            
            # Filter expression file
            real_focal_samples = self.samp_to_class['Sample_name'].tolist()
            self.df_expression = self.df_expression[real_focal_samples]
            
            cap_gene = [str(g).upper() for g in self.df_expression.index.tolist()] # cap the genes
            self.df_expression['Name'] = cap_gene                                  # create a new column
            self.df_expression = self.df_expression[['Name'] + real_focal_samples] # put the 'Name' column at front
            self.df_expression.index = range(0,len(df_expression))                 # number the rows
            
            print "\nFiltered expression file: " + str(self.df_expression.shape)
            display(self.df_expression.head())
            
            self.cls_list = self.samp_to_class[class_name].tolist()
            if check_list(self.cls_list) == False:
                print '\nWarning: Your class column contains only one value. It should contain two. GSEA may not work under these circumstances.'
            
    def call_gsea(self,
              method = 'log2_ratio_of_classes',
              processes = 4,    ## 1 is default
              format = 'png',
              permutation_num = 100, # reduce number to speed up test
              weighted_score_type = 1,  # default: 1
             ):
    
        print "This may take a few minutes."
        self.gs_res = gp.gsea(data = self.df_expression, 
                              gene_sets = self.gmt_file,
                              cls = self.cls_list, 
                              permutation_num = permutation_num,
                              weighted_score_type = weighted_score_type,
                              outdir = self.output_dir,
                              method = method,
                              processes = processes,    ## 1 is default
                              format = format)

        #access the dataframe results throught res2d attribute
        return self.gs_res.res2d.head()
    
    def plot_gsea(self, style_content = 'ggplot', top = 20, y = 'fdr', x = 'Term', fontsize = 8):
        gsea_results = self.gs_res.res2d
        with plt.style.context(style_content):
            gsea_results = gsea_results.reset_index()
            gsea_results.head(top).plot.barh(y = y, x = x, fontsize = fontsize)
    
    def check_list(curr_list):
        first = curr_list[0]
        for obj in curr_list:
            if obj != first:
                return True
        return False