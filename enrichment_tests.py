'''
File name: enrichment_tests.py
Authors: Jimmy Shi, Partha Rao
Date created: 6/28/2017
Date last modified:
Python Version: 2.7
'''

import math
import multiprocessing
import sys
import numpy as np
import time

from scipy import stats
from collections import defaultdict
from flib.core.gmt import GMT
from utilities.enrichment_output_writer import OUT
from utilities.mat import MAT

score_arr = []


# general class for setting up and running an enrichment test
class EnrichmentTest:
    def __init__(self):
        self.test_name=args.test_name
        self.anno_list = GMT(args.annotation_list)
        self.expr_list = MAT(args.expr_list)
        self.expr_cluster = args.cluster_number
        self.permutations = args.permutations
        self.alpha = args.rate
        self.output = args.output
        self.weight = args.weight
        self.cpu_count=args.cpu
        self.precision=args.precision
        self.console=args.console
        self.significant=args.significant

    def run(self):
        '''
        runs an enrichment test, generates an html table

        :return: Nothing, printout() will only write to the specified output file and to the console if specified
        '''

        start_time=time.time()
        rankings = self.switch(self.test_name)
        # passes output to the printer class
        print "Execution time: "+str(time.time()-start_time)
        if self.test_name == "gsea":
            #html table
            OUT(rankings[0], rankings[1], self.output).html_table_GSEA(self.significant,self.precision)
            #output file and console
            return OUT(rankings[0], rankings[1], self.output).printout_GSEA(self.console, self.significant,self.precision)
        #html table
        OUT(rankings[0], rankings[1], self.output).html_table(self.significant,self.precision)
        #output file and console
        return OUT(rankings[0], rankings[1], self.output).printout(self.console, self.significant,self.precision)

    def switch(self, test_name):
        '''
        used by the run method to match the input with the correct test

        :param str test_name: Name of the test to be run
        :return: An array of 2 items, the gene rankings and the filtered FDR rankings
        '''
        if test_name=="wilcoxon":
            return wilcoxon(self.expr_list, self.expr_cluster, self.anno_list, self.alpha,self.cpu_count)
        elif test_name == "page":
            return page(self.expr_list, self.expr_cluster, self.anno_list, self.alpha,self.cpu_count)
        elif test_name == "gsea":
            return gsea(self.expr_list, self.expr_cluster, self.anno_list, self.permutations, self.alpha, self.weight,self.cpu_count)

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


# each enrichment test will require an array of InputItems
# multiprocessing will map each item in the array to a worker
class EnrichmentInputItem:
    def __init__(self, anno_id, anno_list, expr_cluster, expr_list, permutations=None, weight=None, gene_mean=None, gene_sd=None):
        self.anno_id = anno_id
        self.anno_list = anno_list
        self.expr_cluster = expr_cluster
        self.expr_list = expr_list
        self.permutations = permutations
        self.gene_mean = gene_mean
        self.gene_sd = gene_sd
        self.weight=weight

def generate_inputs(anno, expr_cluster, expr_list, permutations=None, weight=None):
    '''
    generates a list of input items to be mapped to several processors

    :param GMT anno: Annotation set object
    :param MAT expr_list: Object containing list of expression values
    :param int expr_cluster: The expression cluster of interest
    :param int permutations: Number of permutations for GSEA, disregarded if test is not GSEA
    :param float weight: The weight of the GSEA test, disregarded if test is not GSEA
    :return: An array of input items
    '''

    input_items = []
    for anno_id in anno.genesets:
        if permutations != None:
            input_items.append(EnrichmentInputItem(anno_id, anno.genesets[anno_id], expr_cluster, expr_list, permutations, weight))
        else:
            input_items.append(EnrichmentInputItem(anno_id, anno.genesets[anno_id], expr_cluster, expr_list))
    return input_items

def multiprocess(input_arr, method, cpu_count):
    '''
    executes a multiprocess which divides the work for each annotation set among processors

    :param list input_arr: The array of input items
    :param func method: The enrichment test to be run
    :param int cpu_count: Number of cores to be used
    :return: The output of a method process --> an EnrichmentResult object
    '''

    p = multiprocessing.Pool(processes=cpu_count)
    results = p.map(method, input_arr)
    p.close()
    p.join()
    return results

def gsea(expr_list, expr_cluster, anno_list, permutations, alpha, weight, cpu_count):
    '''
    the gsea test

    :param MAT expr_list: Object containing list of expression values
    :param int expr_cluster: The expression cluster of interest
    :param GMT anno: Annotation set object
    :param int permutations: Number of permutations for GSEA, disregarded if test is not GSEA
    :param float alpha: The desired significance level
    :param float weight: The weight of the GSEA test, disregarded if test is not GSEA
    :param int cpu_count: Specifies the number of cores to be used
    :return: An array of 2 items, the gene rankings and the filtered FDR rankings
    '''

    expr_list.sort(expr_cluster)
    # build multiprocessing inputs
    input_arr = generate_inputs(anno_list, expr_cluster, expr_list, permutations, weight)

    gene_rankings = multiprocess(input_arr, gsea_process, cpu_count)
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line.p_value))
    return [gene_rankings, significance_filter(gene_rankings, alpha)]

def gsea_process(input_item):
    '''
    gsea sub-method for multiprocessing

    :param list input_item: Object containing input fields
    :return: EnrichmentResult object
    '''

    es = enrichment_score(input_item.anno_list, input_item.expr_cluster, input_item.expr_list, input_item.weight)
    es_arr = es_distr(input_item.expr_list, input_item.expr_cluster, input_item.anno_list, input_item.permutations)
    nes = normalize_score(es, es_arr)
    nes_arr = normalize_array(es_arr)
    n_p = n_p_value(es, es_arr)
    FDR = n_p_value(nes, nes_arr)

    return EnrichmentResult(input_item.expr_cluster, len(input_item.expr_list.dict), input_item.anno_id,
                            len(input_item.anno_list), n_p, FDR, es, nes)

def enrichment_score(anno_set, expr_cluster, expr_list, weight):
    '''
    calculates enrichment score based of the max ES of a linear traversal path

    :param set anno_set: Set of annotation set genes
    :param int expr_cluster: The expression cluster of interest
    :param MAT expr_list: Object containing expression list values
    :param int weight: GSEA weighting option
    :return: The enrichment score
    '''

    anno_map = dict.fromkeys(anno_set, 0)
    anno_map = defaultdict(int, anno_map)
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


def es_distr(expr_list, expr_cluster, anno_list, permutations):
    '''
    generates a distribution of enrichment scores through permutations

    :param MAT expr_list: Object containing list of expression values
    :param int expr_cluster: The expression cluster of interest
    :param GMT anno: Annotation set object
    :param int permutations: Number of permutations for GSEA, disregarded if test is not GSEA
    :return: An array of enrichment scores
    '''

    es_arr = []
    rankings_map = expr_list.dict
    for i in range(0, permutations):
        permuted_arr = np.random.permutation(list(rankings_map))[0:len(anno_list)]
        es_arr.append(enrichment_score(permuted_arr, expr_cluster, expr_list, 1))
    return es_arr



def normalize_score(es, es_arr):
    '''
    normalizes the enrichment score to account for different gene set sizes

    :param int es: Enrichment score
    :param list es_arr: List of enrichment scores obtained through permutations
    :return: The normalized enrichment score
    '''

    total = 0
    count = 0
    for e in es_arr:
        if (e > 0 and es > 0) or (e < 0 and es < 0):
            total += e
            count += 1
    return es / (total / count)



def normalize_array(es_arr):
    '''
    normalizes the array to account for different gene set sizes

    :param list es_arr: List of enrichment scores obtained through permutations
    :return: The normalized array
    '''

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
    '''
    calculates the nominal p value from the enrichment score distribution

    :param int es: Enrichment score
    :param list es_arr: List of enrichment scores obtained through permutations
    :return: The nominal p value
    '''

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


def wilcoxon(expr_list, expr_cluster, anno_list, alpha, cpu_count):
    '''
    the wilcoxon test

    :param MAT expr_list: Object containing list of expression values
    :param int expr_cluster: The expression cluster of interest
    :param GMT anno: Annotation set object
    :param float alpha: The desired significance level
    :param int cpu_count: Specifies the number of cores to be used
    :return: An array of 2 items, the gene rankings and the filtered FDR rankings
    '''

    global score_arr
    score_arr = []
    # build multiprocessing inputs
    input_arr = generate_inputs(anno_list, expr_cluster, expr_list)
    gene_rankings = multiprocess(input_arr, wilcoxon_process, cpu_count)

    rankings_final = benjamini_hochberg(gene_rankings)
    return [rankings_final, significance_filter(rankings_final, alpha)]

def wilcoxon_process(input_item):
    '''
    wilcoxon sub-method for multiprocessing

    :param list input_item: Object containing input fields
    :return: EnrichmentResult object
    '''

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

    return EnrichmentResult(input_item.expr_cluster, len(input_item.expr_list.dict), input_item.anno_id,
                            len(input_item.anno_list), p_value, 0)

def page(expr_list, expr_cluster, anno_list, alpha, cpu_count):
    '''
    the page test

    :param MAT expr_list: Object containing list of expression values
    :param int expr_cluster: The expression cluster of interest
    :param GMT anno: Annotation set object
    :param float alpha: The desired significance level
    :param int cpu_count: Specifies the number of cores to be used
    :return: An array of 2 items, the gene rankings and the filtered FDR rankings
    '''
    score_arr = []

    # calculate value related to the entire cluster
    for gene in expr_list.dict:
        score_arr.append(list(expr_list.dict[gene])[expr_cluster])
    score_arr = np.array(score_arr).astype(np.float)
    gene_mean = np.mean(score_arr)
    gene_sd = np.std(score_arr)
    input_arr = generate_inputs(anno_list, expr_cluster, expr_list)
    input_arr_copy = list(input_arr)

    for i, input_item in enumerate(input_arr_copy):
        # the 0 is included in the input array because of the permutation and weight variables preceding gene_sd and gene_mean in the constructor
        input_arr_copy[i] = EnrichmentInputItem(input_item.anno_id, input_item.anno_list, input_item.expr_cluster,
                                                input_item.expr_list,
                                                0, 0, gene_mean, gene_sd)

    gene_rankings = multiprocess(input_arr_copy, page_process, cpu_count)

    rankings_final = benjamini_hochberg(gene_rankings)
    return [rankings_final, significance_filter(rankings_final, alpha)]

def page_process(input_item):
    '''
    page sub-method for multiprocessing

    :param list input_item: Object containing input fields
    :return: EnrichmentResult object
    '''

    geneset_size = len(input_item.anno_list)
    score_arr = []
    anno_map=dict.fromkeys(input_item.anno_list,0)
    anno_map=defaultdict(int,anno_map)
    # for each gene set, calculate values
    for id in input_item.expr_list.dict:
        row = list(input_item.expr_list.dict[id])
        if id in anno_map:
            score_arr.append(row[input_item.expr_cluster])
    score_arr = np.array(score_arr).astype(np.float)
    geneset_mean = np.mean(score_arr)
    z_score = (input_item.gene_mean - geneset_mean) * math.sqrt(geneset_size) / input_item.gene_sd

    p_value = stats.norm.sf(abs(z_score))
    return EnrichmentResult(input_item.expr_cluster, len(input_item.expr_list.dict), input_item.anno_id,
                            len(input_item.anno_list), p_value, 0)


def benjamini_hochberg(gene_rankings):
    '''
    benjamini hochberg FDR correction for multiple hypothesis testing

    :param list gene_rankings: Object containing EnrichmentResult objects
    :return: A p value adjusted array of EnrichmentResult objects
    '''

    output = []
    # sort rankings before applyign BH correction
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line.p_value))
    prev_bh_value = 0
    for i, E_Result in enumerate(gene_rankings):

        # standard BH correction based on ranking
        bh_value = E_Result.p_value * len(gene_rankings) / float(i + 1)
        bh_value = min(bh_value, 1)
        # ensures that less significant items do not have a higher FDR than more significant items
        if bh_value < prev_bh_value:
            output[i - 1].FDR = bh_value
        E_Result.FDR = bh_value
        output.append(E_Result)
        prev_bh_value = bh_value
    return output

def significance_filter(gene_rankings, alpha):
    '''
    filters out significant items under the alpha level

    :param list gene_rankings: Object containing EnrichmentResult objects
    :param float alpha: The desired significance level
    :return: A filtered array of EnrichmentResult objects
    '''

    significant_values = []
    for i, E_Result in enumerate(gene_rankings):
        if E_Result.FDR <= alpha:
            significant_values.append(E_Result)

    return significant_values


if __name__ == '__main__':
    from argparse import ArgumentParser

    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, version="%prog dev-unreleased")

    parser.add_argument(
        "-n",
        "--name of test to run",
        dest="test_name",
        help="options --> (gsea, page, wilcoxon)",
        metavar="STRING",
        required=True
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
        help="an integer for the number of GSEA gene set permutations (OPTIONAL) (DEFAULT 1000)",
        metavar="INTEGER",
        type=int,
        default=1000
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
        "-w",
        "--GSEA weight",
        dest="weight",
        help="the weighting amount for GSEA (OPTIONAL) (DEFAULT 1)",
        metavar="FLOAT",
        type=float,
        default=1
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
        "--the alpha level for the test",
        dest="rate",
        help="a decimal for the alpha rate (DEFAULT 0.05)",
        metavar="FLOAT",
        type=float,
        default=0.05)
    parser.add_argument(
        "-i",
        "--the number of cores to be used",
        dest="cpu",
        help="an integer for the amount of cores (DEFAULT 1)",
        metavar="INTEGER",
        type=int,
        default=1)
    parser.add_argument(
        "-d",
        "--the decimal precision",
        dest="precision",
        help="an integer for the decimal precision (DEFAULT Full Precision)",
        metavar="INTEGER",
        type=int,
        default=-1)
    parser.add_argument(
        "-t",
        "--print to console",
        dest="console",
        help="a boolean for a print to console option (DEFAULT True)",
        metavar="BOOLEAN",
        type=bool,
        default=True)
    parser.add_argument(
        "-s",
        "--significant only",
        dest="significant",
        help="a boolean for significant only results option (DEFAULT False)",
        metavar="BOOLEAN",
        type=bool,
        default=False)

    args = parser.parse_args()

    # perform a gsea test with a console printout and and including all values
    test = EnrichmentTest()
    test.run()