import math
import multiprocessing
import sys
import time

import numpy as np
from flib.core.gmt import GMT
from scipy import stats

from enrichment_output_writer import OUT
from utilities.mat import MAT
from collections import defaultdict

args=None
score_arr=[]
# builds general multiprocessing inputs for all statistical methods
def generate_inputs(anno, cluster, rankings, permutations=None):

    input_arr = []
    for go_id in anno.genesets:
        if permutations!=None:
            input_arr.append([go_id, anno.genesets[go_id], cluster, rankings,permutations])
        else:
            input_arr.append([go_id, anno.genesets[go_id], cluster, rankings])
    return input_arr


#starts multiprocessing the specified method
def multiprocess(map_arr, method):
    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = p.map(method, map_arr)
    p.close()
    p.join()
    return results


# GSEA multiprocessing function
def gsea_process(m_arr):
    #t1 = time.time()
    es = enrichment_score(m_arr[1], m_arr[2], m_arr[3], 1)
    es_arr = es_distr(m_arr[3], m_arr[2], m_arr[1],m_arr[4])
    nes = normalize_score(es, es_arr)
    nes_arr = normalize_array(es_arr)
    n_p = n_p_value(es, es_arr)
    FDR = n_p_value(nes, nes_arr)
    #print time.time()-t1
    return [m_arr[2], m_arr[0], n_p, FDR, es, nes]


# calculates enrichment score based of the max ES of a linear traversal
def enrichment_score(anno, cluster, rankings, weight):

    anno_map = defaultdict()
    for i in anno:
        anno_map[i]=0
    rankings_map=rankings.dict

    Nhint = len(anno)
    N = len(rankings.dict)
    set_score_sum = 0

    for id in anno:
        if id in rankings_map and len(list(rankings_map[id])) != 0:
            set_score_sum += (float(list(rankings_map[id])[cluster])) ** weight

    max_ES = -sys.maxint
    min_ES = sys.maxint
    net_sum = 0
    sum_arr = []


    for id in rankings.ordered_dict:

        if id in anno_map and len(list(rankings_map[id])) != 0:
            net_sum += (float(list(rankings_map[id])[cluster])) ** weight / set_score_sum
        elif len(list(rankings_map[id])) != 0:
            net_sum += -1.0 / (N - Nhint)

        if net_sum > max_ES:
            max_ES = net_sum
        if net_sum < min_ES:
            min_ES = net_sum
        sum_arr.append(net_sum)
    return max(max_ES,abs(min_ES))


# creates a distribution of enrichment scores from randomization
def es_distr(rankings, cluster, anno,permutations):
    es_arr = []
    rankings_map = rankings.dict
    for i in range(0, permutations):
        permuted_arr = np.random.permutation(list(rankings_map))[0:len(anno)]
        es_arr.append(enrichment_score(permuted_arr, cluster, rankings, 1))
    return es_arr


# gsea main method to call
def gsea(rankings, anno, cluster,permutations, alpha):

    rankings.sort(cluster)
    gene_rankings = []
    input_arr = generate_inputs(anno, cluster, rankings, permutations)
    items = multiprocess(input_arr, gsea_process)

    for i in items:
        gene_rankings.append([i[0], i[1], i[2], i[3], i[4],i[5]])
    # prints out the rankings and significant values
    return [gene_rankings,significance_filter(gene_rankings,alpha)]


def normalize_score(es, es_arr):
    total = 0
    count = 0
    for e in es_arr:
        if (e > 0 and es > 0) or (e < 0 and es < 0):
            total += e
            count += 1
    return es / (total / count)


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
def wilcoxon(rankings, anno, cluster,alpha):

    global score_arr
    score_arr = []

    input_arr=generate_inputs(anno, cluster, rankings)
    gene_rankings= multiprocess(input_arr, wilcoxon_process)

    rankings_final=benjamini_hochberg(gene_rankings)
    return [rankings_final,significance_filter(rankings_final, alpha)]

def wilcoxon_process(m_arr):

    for gene in m_arr[1]:
        row_arr = []
        if gene in m_arr[3].dict:
            row_arr = list(m_arr[3].dict[gene])
        if len(row_arr) != 0:
            score_arr.append(float(row_arr[m_arr[2]]))
        else:
            break
    total_score_list = []
    for gene in m_arr[3].dict:
        row_arr = list(m_arr[3].dict[gene])
        if len(row_arr) != 0 and row_arr[m_arr[2]]:
            total_score_list.append(row_arr[m_arr[2]])
    p_value = stats.ranksums(score_arr, total_score_list)[1]

    return [m_arr[2], m_arr[0], p_value]

# parametric analysis gene enrichment test, compares an input list of genesets versus scores between two experimental groups
def page(rankings, anno, cluster, alpha):
    score_arr = []
    # calculate value related to the entire cluster

    for i in rankings.dict:
        score_arr.append(list(rankings.dict[i])[cluster])
    score_arr = np.array(score_arr).astype(np.float)
    gene_mean = np.mean(score_arr)
    gene_sd = np.std(score_arr)
    input_arr = generate_inputs(anno, cluster, rankings)
    input_arr_copy = list(input_arr)

    for i, row in enumerate(input_arr_copy):
        input_arr_copy[i] = [row[0], row[1], row[2], row[3],
                             gene_mean, gene_sd]

    gene_rankings = multiprocess(input_arr_copy, page_process)

    rankings_final = benjamini_hochberg(gene_rankings)
    return [rankings_final,significance_filter(rankings_final, alpha)]

# page multiprocess method
def page_process(m_arr):

    geneset_size = len(m_arr[1])
    score_arr = []
    # for each gene set, calculate values

    for id in m_arr[3].dict:
        row = list(m_arr[3].dict[id])
        if id in m_arr[1]:
            score_arr.append(row[m_arr[2]])

    score_arr = np.array(score_arr).astype(np.float)
    geneset_mean = np.mean(score_arr)
    z_score = (m_arr[4] - geneset_mean) * math.sqrt(geneset_size) / m_arr[5]

    p_value = stats.norm.sf(abs(z_score))
    return [m_arr[2], m_arr[0], p_value]

# wrapper function to call enrichment tests
def enrichment_test(test_name, print_option):

    anno = GMT(args.annotation_list)
    mat = MAT(args.cluster_list)
    cluster = args.cluster_number
    permutations=args.permutations

    if test_name == "wilcoxon":
        rankings = wilcoxon(mat, anno, cluster,args.rate)
    elif test_name == "page":
        rankings = page(mat, anno, cluster,args.rate)
    elif test_name == "gsea":
        rankings = gsea(mat, anno, cluster, permutations,args.rate)

    # prints out the rankings and significant values

    if test_name == "gsea":
        return OUT(rankings[0],rankings[1], args.output,  mat, anno).printout_GSEA(print_option,False)
    else:
        return OUT(rankings[0],rankings[1], args.output, mat, anno).printout_E(print_option,False)

# FDR correction for multiple hypothesis testing
def benjamini_hochberg(gene_rankings):

    output=[]
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line[2]))
    prev_bh_value = 0
    for i, row in enumerate(gene_rankings):

        bh_value = row[2] * len(gene_rankings) / float(i + 1)
        bh_value = min(bh_value, 1)
        # to preserve monotonicity
        if bh_value < prev_bh_value:
            output[i - 1][3] = bh_value
        temp_arr = []
        for j, grc in enumerate(row):
            temp_arr.append(grc)
        temp_arr.append(bh_value)
        output.append(temp_arr)
        prev_bh_value=bh_value
    return output

#filters out significant items
def significance_filter(gene_rankings,alpha):
    significant_values = []
    for i, row in enumerate(gene_rankings):
        if row[3] < alpha:
            significant_values.append(row)

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
        "-c",
        "--clusters-file",
        dest="cluster_list",
        help="cluster file",
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
        "-l",
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

    # choose fisher_exact, chi_squared, hypergeometric, or binomial
    enrichment_test("gsea", True)
