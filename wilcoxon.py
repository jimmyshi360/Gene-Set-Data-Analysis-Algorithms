from mat import MAT
from flib.core.gmt import GMT
from output_generator import OUT
import scipy.stats
import numpy as np

#wilcoxon rank sum test, compares an input list of genesets versus scores between two experimental groups
def wilcoxon(gmt, mat, output, cluster, false_discovery_rate):
    score_arr = []
    gene_rankings = []

    for gsid in gmt.genesets:
        for gene in gmt.genesets[gsid]:
            row_arr = list(mat.matrix[gene])
            if len(row_arr) != 0:
                score_arr.append(row_arr[cluster])
            total_score_list = []
        for gene in mat.matrix.keys():
            row_arr = list(mat.matrix[gene])
            if len(row_arr) != 0:
                total_score_list.append(row_arr[cluster])
        p_value = scipy.stats.ranksums(score_arr, total_score_list)

        gene_rankings.append([p_value[1], gsid])
        gene_rankings = sorted(gene_rankings, key=lambda line: float(line[0]))

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
    output = open(args.output, "r+")

    wilcoxon(gmt,mat,output, args.cluster_number, args.rate)

