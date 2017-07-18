import math
import multiprocessing
import sys
import time

import numpy as np
from flib.core.gmt import GMT
from scipy import stats

from enrichment_output_writer import OUT
from utilities.mat import MAT
from utilities.parsers import Parsers


# builds general multiprocessing inputs for all statistical methods
def build_inputs(anno, cluster, mat):
    input_arr = []
    for go_id in anno.genesets:
        input_arr.append([go_id, anno.genesets[go_id], cluster, mat.matrix])
    return input_arr


# completes the input array for each sepcific GSEA sample gene set and executes a multiprocess mapping inputs to all GO gene sets
def complete_G_inputs(map_arr, method):
    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = p.map(method, map_arr)
    p.close()
    p.join()
    return results


# GSEA multiprocessing function
def gsea_process(m_arr):
    #t1=time.time()
    es = enrichment_score(m_arr[1], m_arr[2], m_arr[3], 1)
    permuted_arr = permute(m_arr[3], m_arr[1])
    es_arr = es_distr(m_arr[3], m_arr[2], permuted_arr)
    nes = normalize_score(es, es_arr)
    nes_arr = normalize_array(es_arr)
    n_p = n_p_value(es, es_arr)
    FDR = n_p_value(nes, nes_arr)
    #print time.time()-t1
    return [m_arr[2], m_arr[0], n_p, FDR, es, nes]


# calculates enrichment score based of the max ES of a linear traversal
def enrichment_score(anno, cluster, mat, weight):
    Nhint = len(anno)
    N = len(mat)
    set_score_sum = 0

    for id in anno:
        if id in mat and len(list(mat[id])) != 0:
            set_score_sum += (float(list(mat[id])[cluster])) ** weight

    max_ES = -sys.maxint
    min_ES = sys.maxint
    net_sum = 0
    sum_arr = []

    for id in mat:
        if id in anno and len(list(mat[id])) != 0:
            net_sum += (float(list(mat[id])[cluster])) ** weight / set_score_sum
        else:
            net_sum += -1.0 / (N - Nhint)

        if net_sum > max_ES:
            max_ES = net_sum
        if net_sum < min_ES:
            min_ES = net_sum
        sum_arr.append(net_sum)

    return max_ES


# gene set permutation
def permute(mat, anno):
    permuted_arr = []
    for i in range(0, 1000):
        anno_permute = np.random.permutation(list(mat))
        permuted_arr.append(anno_permute[0:len(anno)])
    return permuted_arr


# creates a distribution of enrichment scores from randomization
def es_distr(mat, cluster, permuted_arr):
    es_arr = []
    for i in range(0, len(permuted_arr)):
        es_arr.append(enrichment_score(permuted_arr[i], cluster, mat, 1))
    return es_arr


# gsea main method to call
def gsea(mat, anno, cluster):
    gene_rankings = []
    input_arr = build_inputs(anno, cluster, mat)
    items = complete_G_inputs(input_arr, gsea_process)
    for i in items:
        gene_rankings.append([i[0], i[1], i[2], i[3], i[4]])
    # prints out the rankings and significant values
    return gene_rankings


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
# UNFINISHED
def wilcoxon():
    print("\nWILCOXON")
    p = Parsers("-a -c -o -l -r")
    anno = GMT(p.args.annotation_list)
    mat = MAT(p.args.cluster_list)
    cluster = p.args.cluster_number

    score_arr = []
    gene_rankings = []

    for go_id in anno.genesets:

        for gene in gmt.genesets[gsid]:
            row_arr = list(mat.matrix[gene])
            if len(row_arr) != 0:
                score_arr.append(row_arr[cluster])
        total_score_list = []
        for gene in mat.matrix.keys():
            row_arr = list(mat.matrix[gene])
            if len(row_arr) != 0:
                total_score_list.append(row_arr[cluster])
        p_value = stats.ranksums(score_arr, total_score_list)
        gene_rankings.append([p_value[1], gsid])

    # prints out the rankings and significant values
    return OUT(gene_rankings, p.args.output, p.args.rate).printout_GSEA()


# parametric analysis gene enrichment test, compares an input list of genesets versus scores between two experimental groups
def page(mat, anno, cluster):
    score_arr = []
    # calculate value related to the entire cluster
    for i in mat.matrix:
        score_arr.append(list(mat.matrix[i])[cluster])
    score_arr = np.array(score_arr).astype(np.float)
    gene_mean = np.mean(score_arr)
    gene_sd = np.std(score_arr)
    input_arr = build_inputs(anno, cluster, mat)
    input_arr_copy = list(input_arr)

    for i in range(0, len(input_arr)):
        input_arr_copy[i] = [input_arr_copy[i][0], input_arr_copy[i][1], input_arr_copy[i][2], input_arr_copy[i][3],
                             gene_mean, gene_sd]

    gene_rankings = complete_G_inputs(input_arr_copy, page_process)

    return gene_rankings

# page multiprocess method
def page_process(m_arr):
    geneset_size = len(m_arr[1])
    score_arr = []
    # for each gene set, calculate values
    for id in m_arr[3]:
        row = list(m_arr[3][id])
        if id in m_arr[1]:
            score_arr.append(row[m_arr[2]])

    score_arr = np.array(score_arr).astype(np.float)
    geneset_mean = np.mean(score_arr)
    z_score = (m_arr[4] - geneset_mean) * math.sqrt(geneset_size) / m_arr[5]

    p_value = stats.norm.sf(abs(z_score))

    return [m_arr[2], m_arr[0], p_value]

# wrapper function to call enrichmenet tests
def enrichment_test(test_name, print_option, mat=None, anno=None, rate=None, output=None, cluster=None):
    use_parsers = False

    # if there is no sample input, then parsers are used, sets everything for use with parsers
    if mat == None:
        p = Parsers("-a -c -o -l -r")
        anno = GMT(p.args.annotation_list)
        mat = MAT(p.args.cluster_list)
        cluster = p.args.cluster_number
        use_parsers = True
    else:
        mat = MAT(mat)
        anno = GMT(anno)

    if test_name == "wilcoxon":
        gene_rankings = wilcoxon(mat, anno, cluster)
    elif test_name == "page":
        gene_rankings = page(mat, anno, cluster)
    elif test_name == "gsea":
        gene_rankings = gsea(mat, anno, cluster)

    # prints out the rankings and significant values
    if use_parsers:
        if test_name == "gsea":
            return OUT(gene_rankings, p.args.output, p.args.rate, mat, anno).printout_GSEA(print_option)
        else:
            return OUT(gene_rankings, p.args.output, p.args.rate, mat, anno).printout_E(print_option)

    else:
        if test_name == "gsea":
            return OUT(gene_rankings, output, rate, mat, anno).printout_GSEA(print_option)
        else:
            return OUT(gene_rankings, p.args.output, p.args.rate, mat, anno).printout_E(print_option)


if __name__ == '__main__':
    # choose fisher_exact, chi_squared, hypergeometric, or binomial
    enrichment_test("gsea", True)
