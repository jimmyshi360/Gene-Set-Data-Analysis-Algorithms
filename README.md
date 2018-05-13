# Enrichment tests for genomic data analysis

This repository contains a collection of statistical tests for the enrichment analysis of gene data sets.  This project is fully documented and supported on Python 2.7 compilers. The main methods are located under "elib/core".

The end goal is to connect these methods with an online open source tool for researchers to use. Seven major methods are complete, supporting Cloud computing with room for optimization.

Multiprocessing is supported on local machines with the workload for each method split among "x" number of processors. Currently, annotation gene sets are divided evenly if possible among a multithreaded core.

Link to research presentation: https://docs.google.com/presentation/d/1rApoKuPSkaxc1AxJqH6cGmBFa0JcARlAPCm893c_WME/edit#slide=id.g235264b33b_0_160

![alt text](https://image.ibb.co/fD0i1J/flowchart.png)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Requirements are listed in the requirements.txt

### Installing

Project is pip installable via the github link

## Running the tests

Statistical methods included:
```
Overrepresentation analysis tests:
* Fisher exact
* Binomial
* Chi squared
* Hypergeometric

Data inputs:
* Gene annotations (.gmt)
* Gene set sample (.gmt)
* Optional background genes (.txt)

Enrichment analysis tests:
* Gene Set Enrichment Analysis (GSEA)
* Parametric Analysis of Geneset Expression (PAGE)
* Wilcoxon ranksum test

Data inputs:
* Gene annotations (.gmt)
* Expression list (.mat)
```
## Built With

* [Scipy](https://www.scipy.org/) - Scientific python library
* [Numpy](http://www.numpy.org/) - Mathematical python library

## Versioning

```
Python 2.7
Project version 2.0
```

## Authors

* **Mingze Shi - Lead programmer** 
* **Partha Rao - Assistant programmer for Wilcoxon and Binomial tests** 
See also the list of [contributors](https://github.com/SpecOps167/enrichments/graphs/contributors) who participated in this project.

## License

This project is licensed under the GNU General Public License License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Special thanks to the Troyanskyaya enrichments group for their valuable help with this project

## Setup
```
# setup a virtualenv
virtualenv env
source env/bin/activate

# install dependencies
pip install -r requirements.txt
```
