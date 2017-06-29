from scipy import stats
from flib.core.gmt import GMT
import numpy
import math

# fisher exact p values calculated from 2x2 gene set contingency tables (hypergeometric test substitute)
def fisher_exact(sample, anno, output, false_discovery_rate, background=None):
    gene_rankings = []
    anno_genes = anno._genes
    background_size = 0 if background == None else len(background._genes)
    # contruct and analyze contingency tables

    #2x2 Contigency Table
    ###############LIST########ALL OBSERVED GENES#########
    ######################################################
    #IN ANNO######_______##########_______________########
    #NOT IN ANNO##_______##########_______________########
    ######################################################

    for gsid in sample._genesets:
        list_anno_overlaps = len(sample._genesets[gsid].intersection(anno_genes))
        only_list = len(sample._genesets[gsid]) - list_anno_overlaps
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

sample = GMT("test_files\\ec.topgenes0.1only.ec2_enrich_overlap.gmt")
anno = GMT("test_files\\gobp_human.closed.gmt")
output = open("test_files\\output.txt", "r+")

#overloaded method, you can choose to input a background list too
fisher_exact(sample,anno, output, 0.05)
