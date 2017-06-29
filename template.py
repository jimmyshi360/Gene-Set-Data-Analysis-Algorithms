from mat import MAT
from flib.core.gmt import GMT
from output_generator import OUT
import scipy.stats
import numpy as np

#wilcoxon rank sum test, compares an input list of genesets versus scores between two experimental groups
def test(gmt, mat, anno, output, cluster, false_discovery_rate, background=None):

    gene_rankings = []
    for gsid in gmt.genesets:
        print("obtain values here")
    # prints out the rankings and significant values
    OUT(gene_rankings, output, false_discovery_rate).printout()

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
        required=True
        )
    parser.add_argument(
        "-a",
        "--annotation-list file",
        dest="anno_list",
        help="annotation list file for comparision",
        metavar=".gmt FILE",
        required=True
    )
    parser.add_argument(
        "-c",
        "--clusters-file",
        dest="cluster_list",
        help="cluster file",
        metavar=".mat FILE",
        required=True
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
        "-l",
        "--cluster number",
        dest="cluster_number",
        help="the cluster to run through (DEFAULT 0)",
        metavar="INTEGER",
        type=int,
        default=0)
    parser.add_argument(
        "-r",
        "--false discovery rate",
        dest="rate",
        help="a decimal for the false discovery rate (DEFAULT 0.05)",
        metavar="FLOAT",
        type=float,
        default=0.05)

    args = parser.parse_args()

    gmt = GMT(args.gene_list)
    mat = MAT(args.cluster_list)
    anno= GMT(args.anno_list)
    output = open(args.output, "r+")

    test(gmt,mat, anno, output, args.cluster_number, args.rate)

