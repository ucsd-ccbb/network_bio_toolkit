# network_bio_toolkit

network_bio_toolkit is toolbox of network biology methods currently including Upstream Regulator Analysis functionality, authored by members of the [UCSD Center for Computational Biology & Bioinformatics](http://compbio.ucsd.edu).

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
* visJS2jupyter is a network visualization tool created by us! Visit [here](https://ucsd-ccbb.github.io/visJS2jupyter/) for information about our package.
* mygene is a package we use to translate between different gene naming conventions. Visit [here](http://mygene.info/) for further information.
* ndex2 allows us to store some of our larger databases so you don't have to download any huge files. Just load the networkx graphs right from their servers! Visit [here](http://ndexbio.org/#/) for further information.
* python-louvain (community) is the package we use to perform clustering analysis. Visit [here](https://github.com/taynaud/python-louvain) for further information.
* gprofiler-official is the package we use to annotate our clusters. This package requires Python 3. Visit [here](https://biit.cs.ut.ee/gprofiler/page.cgi?apis) for further information. 
* gseapy is a package that performs gene set enrichment analysis. Visit [here](http://gseapy.readthedocs.io/en/latest/introduction.html) for further information.
* seaborn is a package that makes pretty histograms!  Visit [here](https://seaborn.pydata.org/) for further information.


### Installing

Unfortunately our package is not available yet on Pypi. For now, download our code directly in order to use our toolkit.

## Package summary


## Package tools

### Upstream Regulator Analysis Workflow

module: Upstream

description:

This package includes functions to help you determine potential upstream regulators for a set of differencially expressed genes that you supply. It was inspired by Ingenuity System's [Ingenuity Upstream Regulator Analysis in IPA®](http://pages.ingenuity.com/rs/ingenuity/images/0812%20upstream_regulator_analysis_whitepaper.pdf). Our create_graph module will help you load all of the networks you need for this analysis into networkx graph representations. Our stat_analysis package will calculate p-value and z-score statistics for each potential upstream regulator to determine each one's relevance to your input data set. stat_analysis also includes functions that display information about significant regulators as pandas dataframes or visJS2jupyter networks.

### Network Analysis Workflow

modules: 
- Heat2 (Python 2 compatible version)
- Heat3 (Python 3 compatible version)

description:
- Integrate set of differentially expressed (DE) genes with a user-specified interactome
- User may specify significance cutoffs for inclusion in network propagation
- Interactome may be loaded from NDEX
- Localization analysis of DE genes
- Network propagation from DE gene list
- Clustering of resulting nearest N network proximal genes to input DE genes
- Layout options to focus on individual clusters
- Annotation of clusters (over-representation analysis, with gprofiler)
- Can export to cytoscape for further formatting
- CSV exported with compiled results
- integration with drugbank??

Note:
Annotation of clusters/over-representation analysis is only available for Heat3.py because gprofiler requires python 3. 

### Gene Set Enrichment Analysis Workflow

module: PrepGsea

description:

## Example Notebooks

These are available under the "notebooks" folder.

### ura_notebooks

- URA_Arthritis: This notebook includes all features of the ura module. This notebook is a case study done on a data set produced by the UCSD CCBB. It focuses on the flow of the funcitons in the URA package, with empahsis on the order you should call each function.
- URA_Basic_Example: This notebook provides details on how to use our URA package for your own analyses. 
- URA_Entrez_example: Example notebook for how to run an analysis with entrez genes (gene symbol by default)
- URA_Mouse_example: Example notebook for how to run an analysis with mouse genes (human by default)
- URA_Random_Test: Randomization testing of our URA module
- URA_huvec_brca: This notebook doubles as a validaiton notebook and a real-world use case example. We compare the results of our URA package to the results from the Ingenuity Upstream Regulator Analysis's accompanying paper: [Causal analysis approaches in Ingenuity Pathway Analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3928520/).
- URA_space: Example analysis of space mice genes.  

### clustering_notebooks

- Heat2_example: Basic example notebook fo running Heat2 on a set of mice genes. Includes all features available for our Python 2 network analysis workflow.
- Heat2_GIANT: Example notebook using ndex2 to load [GIANT brain_top](http://giant.princeton.edu/download/) as a networkx graph.
- Heat2_with_ndex2: Example notebook using ndex2 to load the [STRING mouse protein links database](https://string-db.org/cgi/download.pl?sessionId=pFsHqnIzSfjX&species_text=Mus+musculus) as a networkx graph.
- Heat3_with_annotation: Basic example notebook fo running Heat3 on a set of mice genes. Includes all features available for our Python 3 network analysis workflow.

### gsea_notebooks/PrepGsea_package_examples

- PrepGsea_full_example: Example for using PrepGsea. If this is your first time using PrepGsea, use this notebook as your reference. 
- PrepGsea_from_file_example: Example for realoading an already completed analysis back into jupyter notebooks.


## Authors
**Mikayla Webster** (13webstermj@gmail.com)

**Brin Rosenthal** (sbrosenthal@ucsd.edu)

**Mengyi (Miko) Liu** (mikoliu798@gmail.com)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
