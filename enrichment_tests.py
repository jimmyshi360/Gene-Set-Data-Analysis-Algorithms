import math
import multiprocessing
import sys

import numpy as np
from flib.core.gmt import GMT
from scipy import stats

from enrichment_output_writer import OUT
from utilities.mat import MAT
from collections import defaultdict

score_arr = []


# each stat test will return an array of EnrichmentResults
class EnrichmentResult:
    def __init__(self, expr_cluster, expr_list_ngenes, anno_id, anno_ngenes, p_value, FDR, es=None, nes=None):
        self.expr_cluster = expr_cluster
        self.expr_list_ngenes = expr_list_ngenes
        self.anno_id = anno_id
        self.anno_ngenes = anno_ngenes
        self.p_value = p_value
        self.FDR = FDR
        self.es = es
        self.nes = nes

# each stat test will require an array of InputItems
class InputItem:
    def __init__(self, anno_id, anno_list, expr_cluster, expr_list, permutations=None, gene_mean=None, gene_sd=None):
        self.anno_id = anno_id
        self.anno_list = anno_list
        self.expr_cluster = expr_cluster
        self.expr_list = expr_list
        self.permutations = permutations
        self.gene_mean = gene_mean
        self.gene_sd = gene_sd


# builds general multiprocessing inputs for all statistical methods
def generate_inputs(anno, expr_cluster, expr_list, permutations=None):
    input_items = []
    for anno_id in anno.genesets:
        if permutations != None:
            input_items.append(InputItem(anno_id, anno.genesets[anno_id], expr_cluster, expr_list, permutations))
        else:
            input_items.append(InputItem(anno_id, anno.genesets[anno_id], expr_cluster, expr_list))
    return input_items


# starts multiprocessing the specified method
def multiprocess(map_arr, method):
    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = p.map(method, map_arr)
    p.close()
    p.join()
    return results


# gsea main method to call
def gsea(expr_list, expr_cluster, anno_list,  permutations, alpha):
    expr_list.sort(expr_cluster)
    gene_rankings = []
    input_arr = generate_inputs(anno_list, expr_cluster, expr_list, permutations)
    gene_rankings = multiprocess(input_arr, gsea_process)

    # prints out the rankings and significant values
    return [gene_rankings, significance_filter(gene_rankings, alpha)]


# GSEA multiprocessing function
def gsea_process(input_item):
    es = enrichment_score(input_item.anno_list, input_item.expr_cluster, input_item.expr_list, 1)
    es_arr = es_distr(input_item.expr_list, input_item.expr_cluster, input_item.anno_list, input_item.permutations)
    nes = normalize_score(es, es_arr)
    nes_arr = normalize_array(es_arr)
    n_p = n_p_value(es, es_arr)
    FDR = n_p_value(nes, nes_arr)
    return EnrichmentResult(input_item.expr_cluster, len(input_item.expr_list.dict), input_item.anno_id,len(input_item.anno_list), n_p, FDR, es, nes)


# calculates enrichment score based of the max ES of a linear traversal
def enrichment_score(anno_set, expr_cluster, expr_list, weight):
    anno_map = defaultdict()
    for i in anno_set:
        anno_map[i] = 0
    rankings_map = expr_list.dict

    Nhint = len(anno_set)
    N = len(expr_list.dict)
    set_score_sum = 0

    for id in anno_set:
        if id in rankings_map and len(list(rankings_map[id])) != 0:
            set_score_sum += (float(list(rankings_map[id])[expr_cluster])) ** weight

    max_ES = -sys.maxint
    min_ES = sys.maxint
    net_sum = 0

    for id in expr_list.ordered_dict:

        if id in anno_map and len(list(rankings_map[id])) != 0:
            net_sum += (float(list(rankings_map[id])[expr_cluster])) ** weight / set_score_sum
        elif len(list(rankings_map[id])) != 0:
            net_sum += -1.0 / (N - Nhint)

        if net_sum > max_ES:
            max_ES = net_sum
        if net_sum < min_ES:
            min_ES = net_sum

    return max(max_ES, abs(min_ES))


# creates a distribution of enrichment scores through permutations
def es_distr(expr_list, expr_cluster, anno_list, permutations):
    es_arr = []
    rankings_map = expr_list.dict
    for i in range(0, permutations):
        permuted_arr = np.random.permutation(list(rankings_map))[0:len(anno_list)]
        es_arr.append(enrichment_score(permuted_arr, expr_cluster, expr_list, 1))
    return es_arr

#normalizes the enrichment score to account for different gene set sizes
def normalize_score(es, es_arr):
    total = 0
    count = 0
    for e in es_arr:
        if (e > 0 and es > 0) or (e < 0 and es < 0):
            total += e
            count += 1
    return es / (total / count)

#normalizes the array to account for different gene set sizes
def normalize_array(es_arr):
    null_distr = []
    total = 0
    count = 0
    for es in es_arr:
        for e in es_arr:
            if (e > 0 and es > 0) or (e < 0 and es < 0):
                total += e
                count += 1
        null_distr.append(es / (total / count))
    return null_distr

# obtain nominal p value from the enrichment score distribution
def n_p_value(es, es_arr):
    total = np.sum(es_arr)
    mean = total / len(es_arr)
    if es == mean:
        return 1.0
    elif es >= mean:
        tail = []
        for e in es_arr:
            if e >= es:
                tail.append(e)
        return float(len(tail)) / float(len(es_arr))
    else:
        tail = []
        for e in es_arr:
            if e <= es:
                tail.append(e)

        return float(len(tail)) / float(len(es_arr))


# wilcoxon rank sum test, compares an input list of genesets versus scores between two experimental groups
def wilcoxon(expr_list, expr_cluster, anno_list, alpha):
    global score_arr
    score_arr = []

    input_arr = generate_inputs(anno_list, expr_cluster, expr_list)
    gene_rankings = multiprocess(input_arr, wilcoxon_process)

    rankings_final = benjamini_hochberg(gene_rankings)
    return [rankings_final, significance_filter(rankings_final, alpha)]

# wilcoxon multiprocess method
def wilcoxon_process(input_item):
    for gene in input_item.anno_list:
        row_arr = []
        if gene in input_item.expr_list.dict:
            row_arr = list(input_item.expr_list.dict[gene])
        if len(row_arr) != 0:
            score_arr.append(float(row_arr[input_item.expr_cluster]))
        else:
            break
    total_score_list = []
    for gene in input_item.expr_list.dict:
        row_arr = list(input_item.expr_list.dict[gene])
        if len(row_arr) != 0 and row_arr[input_item.expr_cluster]:
            total_score_list.append(row_arr[input_item.expr_cluster])
    p_value = stats.ranksums(score_arr, total_score_list)[1]

    return EnrichmentResult(input_item.expr_cluster, len(input_item.expr_list.dict), input_item.anno_id, len(input_item.anno_list), p_value, 0)


# parametric analysis gene enrichment test
def page(expr_list, expr_cluster, anno_list, alpha):
    score_arr = []

    # calculate value related to the entire cluster
    for i in expr_list.dict:
        score_arr.append(list(expr_list.dict[i])[ expr_cluster])
    score_arr = np.array(score_arr).astype(np.float)
    gene_mean = np.mean(score_arr)
    gene_sd = np.std(score_arr)
    input_arr = generate_inputs(anno_list,  expr_cluster, expr_list)
    input_arr_copy = list(input_arr)

    for i, input_item in enumerate(input_arr_copy):
        # the 0 is included in the input array because of the permutations variable preceding gene_sd and gene_mean in the constructor
        input_arr_copy[i] = InputItem(input_item.anno_id, input_item.anno_list, input_item.expr_cluster, input_item.expr_list,
                                      0, gene_mean, gene_sd)

    gene_rankings = multiprocess(input_arr_copy, page_process)

    rankings_final = benjamini_hochberg(gene_rankings)
    return [rankings_final, significance_filter(rankings_final, alpha)]


# page multiprocess method
def page_process(input_item):
    geneset_size = len(input_item.anno_list)
    score_arr = []

    # for each gene set, calculate values
    for id in input_item.expr_list.dict:
        row = list(input_item.expr_list.dict[id])
        if id in input_item.anno_list:
            score_arr.append(row[input_item.expr_cluster])

    score_arr = np.array(score_arr).astype(np.float)
    geneset_mean = np.mean(score_arr)
    z_score = (input_item.gene_mean - geneset_mean) * math.sqrt(geneset_size) / input_item.gene_sd

    p_value = stats.norm.sf(abs(z_score))
    return EnrichmentResult(input_item.expr_cluster, len(input_item.expr_list.dict), input_item.anno_id,len(input_item.anno_list), p_value, 0)


# wrapper function to call enrichment tests
def enrichment_test(test_name, print_to_console):
    anno_list = GMT(args.annotation_list)
    expr_list = MAT(args.expr_list)
    expr_cluster = args.cluster_number
    permutations = args.permutations

    if test_name == "wilcoxon":
        rankings = wilcoxon(expr_list,  expr_cluster, anno_list, args.rate)
    elif test_name == "page":
        rankings = page(expr_list,  expr_cluster, anno_list,args.rate)
    elif test_name == "gsea":
        rankings = gsea(expr_list,  expr_cluster, anno_list,  permutations, args.rate)

    # prints out the rankings and significant values

    if test_name == "gsea":
        return OUT(rankings[0], rankings[1], args.output).printout_GSEA(print_to_console, False)
    return OUT(rankings[0], rankings[1], args.output).printout_E(print_to_console, False)


# FDR correction for multiple hypothesis testing
def benjamini_hochberg(gene_rankings):
    output = []
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line.p_value))
    prev_bh_value = 0
    for i, E_Result in enumerate(gene_rankings):

        bh_value = E_Result.p_value * len(gene_rankings) / float(i + 1)
        bh_value = min(bh_value, 1)
        # to preserve monotonicity, ensures that comparatively less significant items do not have a more comparatively significant FDR
        if bh_value < prev_bh_value:
            output[i - 1].FDR = bh_value
        E_Result.FDR = bh_value
        output.append(E_Result)
        prev_bh_value = bh_value
    return output


# filters out significant items
def significance_filter(gene_rankings, alpha):
    significant_values = []
    for i, E_Result in enumerate(gene_rankings):
        if E_Result.FDR < alpha:
            significant_values.append(E_Result)

    return significant_values


if __name__ == '__main__':
    from argparse import ArgumentParser

    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, version="%prog dev-unreleased")

    parser.add_argument(
        "-a",
        "--annotations-file",
        dest="annotation_list",
        help="annotation file",
        metavar=".gmt FILE",
        required=True
    )

    parser.add_argument(
        "-e",
        "--gene experessions-file",
        dest="expr_list",
        help="expressions file",
        metavar=".mat FILE",
        required=True
    )
    parser.add_argument(
        "-p",
        "--permutations",
        dest="permutations",
        help="an integer for the number of GSEA permutations (OPTIONAL)",
        metavar="INTEGER",
        type=int

    )
    parser.add_argument(
        "-c",
        "--cluster number",
        dest="cluster_number",
        help="the cluster to run through (DEFAULT 0)",
        metavar="INTEGER",
        type=int,
        default=0)
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
        "--the alpha level for the test",
        dest="rate",
        help="a decimal for the alpha rate (DEFAULT 0.05)",
        metavar="FLOAT",
        type=float,
        default=0.05)

    args = parser.parse_args()

    # perform a gsea test with a console printout
    enrichment_test("gsea", True)
