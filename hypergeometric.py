from scipy import stats
from flib.core.gmt import GMT
import numpy
import math

# fisher exact p values calculated from 2x2 gene set contingency tables (hypergeometric test substitute)
def fisher_exact_test(sample, anno, output, false_discovery_rate, background=None):
    gene_rankings = []
    anno_genome = anno._genes
    background_size = 0 if background == None else len(background._genes)
    # contruct and analyze contingency tables
    for gsid in sample._genesets:
        list_overlaps = len(sample._genesets[gsid].intersection(anno_genome))
        in_list_not_anno = len(sample._genesets[gsid]) - list_overlaps
        not_list_anno = len(anno_genome) - list_overlaps
        not_in_anno_and_list = background_size - list_overlaps - in_list_not_anno - not_list_anno if background != None else 0
        next_p_value = stats.fisher_exact([[list_overlaps, not_list_anno], [in_list_not_anno, not_in_anno_and_list]])[1]
        gene_rankings.append([next_p_value, gsid])

    # sort array in descending order
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line[0]))

    # grouping significant gene sets and mul
    # tiple hypothesis test correction (hochberg)
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

sample = GMT("C:\\Users\\Jimmy\\Documents\\GENOMICS\\list.txt")
anno = GMT("C:\\Users\\Jimmy\\Documents\\GENOMICS\\anno.txt")
output = open("C:\\Users\\Jimmy\\Documents\\GENOMICS\\output.txt", "r+")

fisher_exact_test(sample,anno, output, 0.05)
