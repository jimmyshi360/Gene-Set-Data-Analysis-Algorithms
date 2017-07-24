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

class EnrichmentResult:
    def __init__(self, cluster, cluster_ngenes, go_id, gs2_ngenes, p_value, FDR,es=None, nes=None):
        self.cluster=cluster
        self.cluster_ngenes=cluster_ngenes
        self.go_id = go_id
        self.gs2_ngenes=gs2_ngenes
        self.p_value = p_value
        self.FDR=FDR
        self.es=es
        self.nes=nes

    def set_FDR(self,new_FDR):
        self.FDR=new_FDR

class InputItem:
    def __init__(self, go_id, gene_list, cluster,rankings, permutations=None, gene_mean=None, gene_sd=None):
        self.go_id = go_id
        self.gene_list=gene_list
        self.cluster=cluster
        self.rankings=rankings
        self.permutations= permutations
        self.gene_mean=gene_mean
        self.gene_sd=gene_sd

# builds general multiprocessing inputs for all statistical methods
def generate_inputs(anno, cluster, rankings, permutations=None):

    input_arr = []
    for go_id in anno.genesets:
        if permutations!=None:
            input_arr.append(InputItem(go_id, anno.genesets[go_id], cluster, rankings,permutations))
        else:
            input_arr.append(InputItem(go_id, anno.genesets[go_id], cluster, rankings))
    return input_arr


#starts multiprocessing the specified method
def multiprocess(map_arr, method):
    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = p.map(method, map_arr)
    p.close()
    p.join()
    return results

# gsea main method to call
def gsea(rankings, anno, cluster,permutations, alpha):

    rankings.sort(cluster)
    gene_rankings = []
    input_arr = generate_inputs(anno, cluster, rankings, permutations)
    gene_rankings = multiprocess(input_arr, gsea_process)

    # prints out the rankings and significant values
    return [gene_rankings,significance_filter(gene_rankings,alpha)]

# GSEA multiprocessing function
def gsea_process(input_item):
    #t1 = time.time()
    es = enrichment_score(input_item.gene_list, input_item.cluster, input_item.rankings, 1)
    es_arr = es_distr(input_item.rankings, input_item.cluster, input_item.gene_list, input_item.permutations)
    nes = normalize_score(es, es_arr)
    nes_arr = normalize_array(es_arr)
    n_p = n_p_value(es, es_arr)
    FDR = n_p_value(nes, nes_arr)
    #print time.time()-t1
    return EnrichmentResult(input_item.cluster, len(input_item.rankings.dict), input_item.go_id, len(input_item.gene_list), n_p, FDR, es, nes)

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

def wilcoxon_process(input_item):

    for gene in input_item.gene_list:
        row_arr = []
        if gene in input_item.rankings.dict:
            row_arr = list(input_item.rankings.dict[gene])
        if len(row_arr) != 0:
            score_arr.append(float(row_arr[input_item.cluster]))
        else:
            break
    total_score_list = []
    for gene in input_item.rankings.dict:
        row_arr = list(input_item.rankings.dict[gene])
        if len(row_arr) != 0 and row_arr[input_item.cluster]:
            total_score_list.append(row_arr[input_item.cluster])
    p_value = stats.ranksums(score_arr, total_score_list)[1]

    return EnrichmentResult(input_item.cluster, len(input_item.rankings.dict), input_item.go_id, len(input_item.gene_list), p_value,0)

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

    for i, input_item in enumerate(input_arr_copy):

        #the 0 is included in the input array because of the permutations variable preceding gene_sd and gene_mean in the constructor
        input_arr_copy[i] = InputItem(input_item.go_id, input_item.gene_list, input_item.cluster, input_item.rankings, 0, gene_mean, gene_sd)

    gene_rankings = multiprocess(input_arr_copy, page_process)

    rankings_final = benjamini_hochberg(gene_rankings)
    return [rankings_final,significance_filter(rankings_final, alpha)]

# page multiprocess method
def page_process(input_item):

    geneset_size = len(input_item.gene_list)
    score_arr = []
    # for each gene set, calculate values

    for id in input_item.rankings.dict:
        row = list(input_item.rankings.dict[id])
        if id in input_item.gene_list:
            score_arr.append(row[input_item.cluster])

    score_arr = np.array(score_arr).astype(np.float)

    geneset_mean = np.mean(score_arr)

    z_score = (input_item.gene_mean - geneset_mean) * math.sqrt(geneset_size) / input_item.gene_sd

    p_value = stats.norm.sf(abs(z_score))
    return EnrichmentResult(input_item.cluster, len(input_item.rankings.dict), input_item.go_id, len(input_item.gene_list),
                            p_value, 0)

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
        return OUT(rankings[0],rankings[1], args.output).printout_GSEA(print_option,False)
    else:
        return OUT(rankings[0],rankings[1], args.output).printout_E(print_option,False)

# FDR correction for multiple hypothesis testing
def benjamini_hochberg(gene_rankings):

    output=[]
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line.p_value))
    prev_bh_value = 0
    for i, E_Result in enumerate(gene_rankings):

        bh_value = E_Result.p_value * len(gene_rankings) / float(i + 1)
        bh_value = min(bh_value, 1)
        # to preserve monotonicity
        if bh_value < prev_bh_value:
            output[i - 1].set_FDR (bh_value)
        E_Result.set_FDR(bh_value)
        output.append(E_Result)
        prev_bh_value=bh_value
    return output

#filters out significant items
def significance_filter(gene_rankings,alpha):
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
