# Enrichments tests for genomic data analysis

One Paragraph of project description goes here

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

Requirements are listed in the requirements.txt

### Installing

Project is pip installable via the github link

## Running the tests

Statistical methods included:

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

## Built With

* [Scipy](https://www.scipy.org/) - Scientific python library
* [Numpy](http://www.numpy.org/) - Mathematical python library

## Versioning

Version 2.0

## Authors

* **Mingze Shi** 
* **Partha Rao** 
See also the list of [contributors](https://github.com/SpecOps167/enrichments/graphs/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

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
