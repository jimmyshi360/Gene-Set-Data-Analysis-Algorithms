from mat import MAT
from flib.core.gmt import GMT
import scipy.stats
import numpy as np

#wilcoxon rank sum test, compares an input list of genesets versus scores between two experimental groups
def wilcoxon(gmt, mat, output, cluster, false_discovery_rate):
    score_arr = []
    gene_rankings = []

    for gsid in gmt._genesets:
        for gene in gmt._genesets[gsid]:
            row_arr = list(mat._matrix[gene])
            if len(row_arr) != 0:
                score_arr.append(row_arr[cluster])
        total = mat._matrix.keys()
        total_score_list = []
        for gene in total:
            row_arr = list(mat._matrix[gene])
            if len(row_arr) != 0:
                total_score_list.append(row_arr[cluster])
        p_value = scipy.stats.ranksums(score_arr, total_score_list)

        gene_rankings.append([p_value[1], gsid])
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

