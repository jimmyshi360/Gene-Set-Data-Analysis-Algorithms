# Gene Set Data Analysis Enrichment Statistical Tests

A collection of statistical tests for the enrichment analysis of gene data sets.  Main methods are located under "elib/core".

The end goal is to connect these methods with an online tool for quick and easy enrichment analysis of gene fold data. Functions are optimized for efficiency leveraging python's multiprocessing library.

Since I like pretty pictures: 

![alt text](https://image.ibb.co/fD0i1J/flowchart.png)

Link to research presentation: https://docs.google.com/presentation/d/1rApoKuPSkaxc1AxJqH6cGmBFa0JcARlAPCm893c_WME/edit#slide=id.p

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Installation

Project is pip installable via the github link

## Example Usage

### Example Usage

```python
  enrichment_test = EnrichmentTest(name, GMT_list, MAT_list, clusters, permutations,  alpha, output, weight, cpu_count);
  enrichment_test.run()
```

Or run it in a CLI with:

```
  python enrichment_tests.py -n name -a GMT_list -e MAT_list -c clusters -p permutations -r alpha -o output -w weight -i cpu_count
```

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
Project version 2.1.0
```

## Authors

* **Jimmy Shi - Lead programmer** 
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
