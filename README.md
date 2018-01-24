# network_bio_toolkit

network_bio_toolkit is toolbox of network biology methods currently including Upstream Regulator Analysis functionality, authored by members of the [UCSD Center for Computational Biology & Bioinformatics](http://compbio.ucsd.edu).

## Getting Started

These instructions will get you a copy of the package up and running on your local machine.

### Prerequisites

You must have Jupyter Notebook already installed. Visit [here](http://jupyter.org/install.html) for more information.

The following packages must be installed before using network_bio_toolkit, though most are installed by default if you installed Jupyter Notebook using Anaconda:

* scipy
* pandas
* numpy
* networkx
* matplotlib
* mygene
* visJS2jupyter

If you do not already have these packages by default:
* To install matplotlib or for further information, visit [here](http://matplotlib.org/users/installing.html).
* To install networkx or for further information, visit [here](https://networkx.github.io/)
 
Packages you most likely will not already have:
* mygene is a package we use to translate between different gene naming conventions. Visit [here](http://mygene.info/) for further information.
* visJS2jupyter is a network visualization package created by the UCSD CCBB. Visit [here](https://ucsd-ccbb.github.io/visJS2jupyter/) for information about our package.
 
To install these packages in your comand prompt, use the following commands:

```
>> pip install matplotlib
>> pip install networkx
>> pip install mygene
>> pip install visJS2jupyter
```


### Installing

Information on importing network_bio_toolkit coming soon.


## Features and Examples

### Upstream Regulator Analysis (ura) package
This package includes functions to help you determine potential upstream regulators for a set of differencially expressed genes that you supply. It was inspired by Ingenuity System's [Ingenuity Upstream Regulator Analysis in IPA®](http://pages.ingenuity.com/rs/ingenuity/images/0812%20upstream_regulator_analysis_whitepaper.pdf). Our create_graph module will help you load all of the networks you need for this analysis into networkx graph representations. Our stat_analysis package will calculate p-value and z-score statistics for each potential upstream regulator to determine each one's relevance to your input data set. stat_analysis also includes functions that display information about significant regulators as pandas dataframes or visJS2jupyter networks.

### Example Notebooks

These are available under the "notebooks" folder

1) **URA_Basic_Example**: This notebook provides details on how to use our URA package for your own 
2) **URA_huvec_brca**: This notebook doubles as a validaiton notebook and a real-world use case example. We compare the results of our URA package to the results from the Ingenuity Upstream Regulator Analysis's accompanying paper: [Causal analysis approaches in Ingenuity Pathway Analysis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3928520/).
3) **URA_Arthritis**: This notebook is a case study done on a data set produced by the UCSD CCBB. It focuses on the flow of the funcitons in the URA package, with empahsis on the order you should call each function.


## Authors
**Mikayla Webster** (m1webste@ucsd.edu)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details