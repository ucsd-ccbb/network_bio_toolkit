# network_bio_toolkit

network_bio_toolkit is toolbox of network biology methods authored by members of the [UCSD Center for Computational Biology & Bioinformatics](http://compbio.ucsd.edu).

## Package Summary

Network biology is quickly becoming an indispensable tool for elucidating dysregulated biological pathways and furthering understanding in disease. However, network analysis is intimidating to many researchers.  We present an integrated network analysis toolkit, network_bio_toolkit, built upon our previous interactive network visualization tool, visJS2jupyter, to provide three common network analysis workflows for RNAseq data, from start to finish, requiring only an expression file, or a list of differentially expressed genes.  The tool is designed for use in Jupyter notebooks, to promote reproducible research.

This package includes a set of start-to-finish functions for three common network analysis pipelines: upstream regulator analysis, network propagation and clustering, and gene set enrichment analysis. The upstream regulator analysis pipeline includes localization, transcription factor enrichment, and transcription factor activation state prediction methods. The network propagation and clustering pipeline includes network propagation, clustering, and annotation functionality. The gene set enrichment analysis includes easy-to-use data filtering methods, that help the user prep their data as input for GSEA’s enrichment calculation function. Each pipeline includes functions to easily load and process a variety of genome networks, such as STRING and GIANT, from NDEXbio, eliminating the need for downloading large and unwieldy network files. Each pipeline also has data visualization functions, so interpreting results from our pipelines is intuitive. Many of the analysis components of each pipeline are modular; the user can pick and choose which functions are relevant to their data and leave out unnecessary steps. We demonstrate the functionality of network_bio_toolkit by analyzing arthritis, breast cancer, HUVEC, and space mice differential expression gene sets using our three pipelines.

## Getting Started

These instructions will get a copy of the package up and running on your local machine.

### Prerequisites

All modules (expect for Heat3.py) were created to work with Python 2. We can't guarentee that all modules will work with Python 3.

You must have the following packages installed in order to use network_bio_toolkit.

The following are common packages, that are most likely already installed on your system:

* pandas
* scipy
* networkx
* numpy
* matplotlib

Please note that our code expects networkx version 1.11. We do not guarentee that all code will work with networkx 2.1 or greater. To install networkx version 1.11, run the following code:

```
>> pip install networkx==1.11
```

Further references:
* For further information on matplotlib, visit [here](http://matplotlib.org/users/installing.html).
* For further information on networkx, visit [here](https://networkx.github.io/).

You must also have the following packages installed:

* visJS2jupyter
* mygene
* ndex2
* community
* gprofiler
* gseapy
* seaborn

To install these packages, please run the following commands, respectively.

```
>> pip install visJS2jupyter
>> pip install mygene
>> pip install ndex2
>> pip install python-louvain
>> pip install gprofiler-official
>> pip install gseapy
>> pip install seaborn
```

Further references:
* **visJS2jupyter** is a network visualization tool created by us! Visit [here](https://ucsd-ccbb.github.io/visJS2jupyter/) for information about our package.
* **mygene** is a package we use to translate between different gene naming conventions. Visit [here](http://mygene.info/) for further information.
* **ndex2** allows us to store some of our larger databases so you don't have to download any huge files. Just load the networkx graphs right from their servers! Visit [here](http://ndexbio.org/#/) for further information.
* **python-louvain** (community) is the package we use to perform clustering analysis. Visit [here](https://github.com/taynaud/python-louvain) for further information.
* **gprofiler-official** is the package we use to annotate our clusters. This package requires Python 3. Visit [here](https://biit.cs.ut.ee/gprofiler/page.cgi?apis) for further information. 
* **gseapy** is a package that performs gene set enrichment analysis. Visit [here](http://gseapy.readthedocs.io/en/latest/introduction.html) for further information.
* **seaborn** is a package that makes pretty histograms!  Visit [here](https://seaborn.pydata.org/) for further information.


### Installing

Either download our package from github or pip install from pypi.

```
pip install network_bio_toolkit
```

## Package tools

### Upstream Regulator Analysis Workflow

**module**: Upstream

**description**: This package includes functions to help you determine potential upstream regulators for a set of differentially expressed genes that you supply. 

**workflow**:

- Integrate set of differentially expressed (DE) genes with a user-specified interactome.  Interactome may be loaded from NDEX (no downloading required)
- Seed genes for network analysis defined by user-specified significance cutoffs in DE gene file.  User may select a number of options here, including whether to use gene symbol, or entrez id, and which columns to use for significance cutoffs and fold change cutoffs. 
- Network localization analysis of DE genes, using one of two methods.  The first method, 'num_edges', calculates the number of edges contained in the subgraph composed of a sub-sampling of the DE genes, and compares to degree-matched, randomly selected genes.  This method is similar to the number of edges p-value used in STRING-DB. The second method, 'LCC', calculates the size of the largest connected component in the DE genes subgraph, compared to degree-matached randomly selected genes.  This method was used in Menche et al 2015 (http://science.sciencemag.org/content/347/6224/1257601).  Sampling and randomization are repeated a user-specified number of times (we recommend > 100).  Network localization allows for determining if the input seed gene set is more connected in the network than would be expected by chance, and thus likely representing an underlying biological mechanism/pathway.  
- Transcription factor (TF) enrichment enrichment calculations, using a hypergeometric test (calculates the log of the p-value)
- TF activation state prediction, represented as a signed z-score. TFs assigned positive z-scores are activating, while negative z-scores are inhibiting. Significance of z-score indicates the confidence of the activation state prediction.
- Visualization methods to pinpoint most significant activating and inhbiting TFs, display distribution of TFs as a histogram, and display individual TF-plus-downstream-target networks using visJS2jupyter.

### Heat Propagation and Clustering Analysis Workflow

**modules**: 
- Heat2 (Python 2 compatible version)
- Heat3 (Python 3 compatible version)

**description**: This package provides tools to conduct an integrated network analysis of a set of differentially expressed genes.  

**workflow**:

- Integrate set of differentially expressed (DE) genes with a user-specified interactome.  Interactome may be loaded from NDEX (no downloading required)
- Seed genes for network analysis defined by user-specified significance cutoffs in DE gene file.  User may select a number of options here, including whether to use gene symbol, or entrez id, and which columns to use for significance cutoffs and fold change cutoffs. 
- Network localization analysis of DE genes, using one of two methods.  The first method, 'num_edges', calculates the number of edges contained in the subgraph composed of a sub-sampling of the DE genes, and compares to degree-matched, randomly selected genes.  This method is similar to the number of edges p-value used in STRING-DB. The second method, 'LCC', calculates the size of the largest connected component in the DE genes subgraph, compared to degree-matached randomly selected genes.  This method was used in Menche et al 2015 (http://science.sciencemag.org/content/347/6224/1257601).  Sampling and randomization are repeated a user-specified number of times (we recommend > 100).  Network localization allows for determining if the input seed gene set is more connected in the network than would be expected by chance, and thus likely representing an underlying biological mechanism/pathway.  
- Network propagation from DE gene list.  Network propagation provides information about the local neighborhood of a set of ‘seed’ genes in network space, and genes that have high network propagation scores are likely to be related to these ‘seed’ genes.  See https://www.nature.com/articles/nrg.2017.38 for a review.
- Clustering of resulting nearest N (we recommend N = 300-500) network proximal genes to input DE genes, forming a 'proximal subnetwork'.  Clustering allows for exploration and annotation of the pathways represented in the proximal subnetwork.  Clustering is conducted using a network modularity maximization algorithm (https://github.com/taynaud/python-louvain, http://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008/meta), which identifies groups of nodes which are more highly connected within the group than between groups.  Annotation of clusters is performed using over-representation analysis, with gprofiler, using as the background set the N genes in the proximal subnetwork.  
- We additionally provide layout options to focus on individual clusters, a table containing the full analysis results, and export options to cytopscape for further formatting.  

**Note**: Annotation of clusters/over-representation analysis is only available for Heat3.py because gprofiler requires python 3. 

### Gene Set Enrichment Analysis Workflow

**module**: PrepGsea

**description**: The gene set enrichment analysis workflow provides an interactive interface to the gseapy tools, consistent with the format used in our other tools. 

**workflow**:

- The user provides an expression file and an experiment metadata file
- The user selects the desired comparison for GSEA from provided interactive widgets. 
- The resulting enriched pathways and GSEA figures are saved in a specified output folder (we also provide tools to load pre-computed results)
- The user may plot the enriched pathways in barcharts and heatmaps for inspection of individually up or down regulated genes.    

## Example Notebooks

These are available under the "notebooks" folder.

### ura_notebooks

- **URA_Arthritis**: This notebook includes all features of the ura module. This notebook is a case study done on a data set produced by the UCSD CCBB. It focuses on the flow of the funcitons in the URA package, with empahsis on the order you should call each function.
- **URA_Basic_Example**: This notebook provides details on how to use our URA package for your own analyses. 
- **URA_Entrez_example**: Example notebook for how to run an analysis with entrez genes (gene symbol by default)
- **URA_Mouse_example**: Example notebook for how to run an analysis with mouse genes (human by default)
- **URA_Random_Test**: Randomization testing of our URA module
- **URA_huvec_brca**: This notebook doubles as a validaiton notebook and a real-world use case example. We compare the results of our URA package to the results from the Ingenuity Upstream Regulator Analysis's accompanying paper: [Causal analysis approaches in Ingenuity Pathway Analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3928520/).
- **URA_space**: This analysis was conducted on the liver genes of NASA lab mice who had experienced space flight. Space mice!

### clustering_notebooks

- **Heat2_example**: Basic example notebook fo running Heat2 on a set of mice genes. Includes all features available for our Python 2 network analysis workflow.
- **Heat2_GIANT**: Example notebook using ndex2 to load [GIANT brain_top](http://giant.princeton.edu/download/) as a networkx graph.
- **Heat2_with_ndex2**: Example notebook using ndex2 to load the [STRING mouse protein links database](https://string-db.org/cgi/download.pl?sessionId=pFsHqnIzSfjX&species_text=Mus+musculus) as a networkx graph.
- **Heat3_with_annotation**: Basic example notebook fo running Heat3 on a set of mice genes. Includes all features available for our Python 3 network analysis workflow.

### gsea_notebooks

- **PrepGsea_basic_example**: Example for using PrepGsea. If this is your first time using PrepGsea, use this notebook as your reference. 
- **PrepGsea_from_file_example**: Example for realoading an already completed analysis back into jupyter notebooks.

### general_notebooks

- Look [here](https://bl.ocks.org/m1webste/raw/150c2f9f1c00209f02cbd0e6c504de60/) for an example of all three pipelines in one notebook. This analysis was conducted on the kidney genes of NASA lab mice who had experienced space flight. The actual notebook that bl.ocks page was derived from is available in our notebooks/general_notebooks folder.


## Misc

### Background Database Sources for Download

Our packages include functions for loading STRING's protein actions and protein links files for mouse and human genomes.
- [STRING Homo Sapiens](https://string-db.org/cgi/download.pl?sessionId=JSIcPdRainGc&species_text=Homo+sapiens)
- [STRING Mus Musculus](https://string-db.org/cgi/download.pl?sessionId=pFsHqnIzSfjX&species_text=Mus+musculus)

### Background Databases from ndex2
Our full set of ndex databases is available [here](http://ndexbio.org/#/user/9f248194-480b-11e8-a935-0ac135e8bacf).

## Authors
- **Mikayla Webster** (13webstermj@gmail.com)
- **Brin Rosenthal** (sbrosenthal@ucsd.edu)
- **Mengyi (Miko) Liu** (mikoliu798@gmail.com)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
