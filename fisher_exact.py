from scipy import stats
from flib.core.gmt import GMT
import numpy
import math


# fisher exact test, compares an input list of genes against an annotation list and all known/observed genes
def fisher_exact(sample, anno, output, false_discovery_rate, background=None):
    gene_rankings = []
    anno_genes = anno._genes
    background_size = 0 if background == None else len(background._genes)
    # contruct and analyze contingency tables

    #2x2 Contigency Table
    ###############LIST########ALL OBSERVED GENES########
    #####################################################
    #IN ANNO######_______##########_______________#######
    #NOT IN ANNO##_______##########_______________#######
    #####################################################

    for gsid in sample.genesets:
        list_anno_overlaps = len(sample.genesets[gsid].intersection(anno_genes))
        only_list = len(sample.genesets[gsid]) - list_anno_overlaps
        anno_only = len(anno_genes) - list_anno_overlaps
        genome_only = background_size - list_anno_overlaps - only_list - anno_only if background != None else 0
        p_value = stats.fisher_exact([[list_anno_overlaps, anno_only], [only_list, genome_only]])[1]
        gene_rankings.append([p_value, gsid])

    # sort array in descending order
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line[0]))

    # grouping significant gene sets and multiple hypothesis test correction (hochberg)
    significant_values = []
    for i in range(0, len(gene_rankings)):
        if gene_rankings[i][0] < float(i) / len(gene_rankings) * false_discovery_rate:
            significant_values.append(gene_rankings[i])
    # print rankings and write to output file, reverses ascending array before sorting
    print("\n\nRANKINGS")
    for set_arr in gene_rankings[::-1]:
        output.write(set_arr[1] + ": " + str(set_arr[0]))
        print(set_arr[1] + ": " + str(set_arr[0]))

    # print all significant gene sets
    print("\n\nSIGNIFICANT VALUES")
    for x in significant_values:
        print(x[1] + " " + str(x[0]))

if __name__ == '__main__':
    from argparse import ArgumentParser

    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, version="%prog dev-unreleased")

    parser.add_argument(
        "-g",
        "--gene-list file",
        dest="gene_list",
        help="gene list file for comparision",
        metavar=".gmt FILE",
        required = True
        )
    parser.add_argument(
        "-a",
        "--annotations-file",
        dest="annotation_list",
        help="annotation file",
        metavar=".gmt FILE",
        required=True
        )
    parser.add_argument(
        "-b",
        "--background-gene file",
        dest="background_list",
        help="background gene list file for comparision (OPTIONAL)",
        metavar=".gmt FILE"
        )
    parser.add_argument(
        "-o",
        "--output-file",
        dest="output",
        help="file to output",
        metavar=".txt FILE",
        required=True
        )
    parser.add_argument(
        "-r",
        "--false discovery rate",
        dest="rate",
        help="a decimal for the false discovery rate (DEFAULT 0.05)",
        metavar="FLOAT",
        type=float,
        default=0.05)

    args = parser.parse_args()
    sample = GMT(args.gene_list)
    anno = GMT(args.annotation_list)
    output = open(args.output,"r+")
    #overloaded method, you can choose to input a background list too
    fisher_exact(sample,anno, output, args.rate)
