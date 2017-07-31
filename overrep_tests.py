'''
File name: overrep_tests.py
Authors: Jimmy Shi, Partha Rao
Date created: 6/28/2017
Date last modified:
Python Version: 2.7
'''

import multiprocessing
import time

from flib.core.gmt import GMT
from scipy import stats
from utilities.background import BACKGROUND
from utilities.overrep_ouput_writer import OUT

'''
2x2 contingency table for over representation tests

______________________________Gene List_______________________Genome________________
In Annotations     |    list_anno_overlaps    |       genome_anno_overlaps         |
------------------------------------------------------------------------------------
Not in Annotations |    list_genome_overlaps  |          genome_only               |
------------------------------------------------------------------------------------

http://bg.upf.edu/course-bco/documents/enrichment.pdf

MULTIPROCESSING:
-Passes in a map_arr storing a list of inputs 
-Used to condense GO term comparison loop

must first call build_inputs() to generate data inputs
call multiprocess() to generate output
'''


# general class for setting up and running an over-representation test

class OverrepTest():
    def __init__(self):
        self.test_name=args.test_name
        self.sample_sets = GMT(args.gene_sets)
        self.anno_list = GMT(args.annotation_list)
        self.background = BACKGROUND([], args.background_list)
        self.alpha = args.rate
        self.output = args.output
        self.cpu_count = args.cpu
        self.precision = args.precision
        self.console=args.console
        self.significant=args.significant

    def run(self):
        '''
        runs an overrep test, generates an HTML table

        :return: Nothing, printout() will only write to the specified output file and to the console if specified
        '''

        start_time = time.time()
        rankings = self.switch(self.test_name)

        print "Execution time: "+str(time.time()-start_time)
        #html table
        OUT(rankings[0], rankings[1], self.output).html_table(self.significant, self.precision)
        #output file and table
        return OUT(rankings[0], rankings[1], self.output).printout(self.console, self.significant, self.precision)

    def switch(self, test_name):
        '''
        used by the run method to match the input with the correct test

        :param str test_name: The name of the test to be executed
        :return: An array of 2 items, the gene rankings and the filtered FDR rankings
        '''

        if test_name == "fisher_exact":
            return fisher_exact(self.sample_sets, self.anno_list, self.alpha, self.background, self.cpu_count)
        elif test_name == "chi_squared":
            return chi_squared(self.sample_sets, self.anno_list, self.alpha, self.background, self.cpu_count)
        elif test_name == "binomial":
            return binomial(self.sample_sets, self.anno_list, self.alpha, self.background, self.cpu_count)
        elif test_name == "hypergeometric":
            return hypergeometric(self.sample_sets, self.anno_list, self.alpha, self.background, self.cpu_count)


# each stat test will return an array of OverrepResults
# contains all useful information to be outputted
class OverrepResult:
    def __init__(self, gsid, sample_set_ngenes, anno_id, anno_ngenes, p_value, overlaps, FDR):
        self.gsid = gsid
        self.sample_set_ngenes = sample_set_ngenes
        self.anno_id = anno_id
        self.anno_ngenes = anno_ngenes
        self.p_value = p_value
        self.overlaps = overlaps
        self.FDR = FDR


# each overrep test will require an array of InputItems
# multiprocessing will map each item in the array to a worker
class OverrepInputItem:
    def __init__(self, anno_id, anno_list, background, gsid, gene_set):
        self.anno_id = anno_id
        self.anno_list = anno_list
        self.background = background
        self.gsid = gsid
        self.gene_set = gene_set


def generate_inputs(anno, background):
    '''
    generates a list of incomplete input items to be mapped to several processors
    :param int gsid: Id of the current gene set
    :param GMT sample: Sample set object
    :param BACKGROUND background: Object containing background gene information
    :return: The output of a method process --> an OverrepResult object
    '''

    input_items = []
    for anno_id in anno.genesets:
        input_items.append([anno_id, anno.genesets[anno_id], background])
    return input_items

def multiprocess(gsid, sample, input_arr, method, cpu_count):
    '''
    completes the input array for each specific sample gene set
    each gene set has its own additional input items to be appended to the array
    executes a multiprocess which divides the work for each annotation set among processors

    :param int gsid: Id of the current gene set
    :param GMT sample: Sample set object
    :param list input_arr: Incomplete input array
    :param func method: The overrep test to be run
    :param int cpu_count: The number of cores to be used
    :return: The output of a method process --> an OverrepResult object
    '''

    input_items = []
    for i, row in enumerate(input_arr):
        anno_id = row[0]
        anno_list = row[1]
        background = row[2]

        # gsid and the gene set are the above mentioned additional items
        next_item = OverrepInputItem(anno_id, anno_list, background, gsid, sample.genesets[gsid])
        input_items.append(next_item)

    p = multiprocessing.Pool(processes=cpu_count)
    results = p.map(method, input_items)
    p.close()
    p.join()

    return results


def gen_table(sample_set, anno_set, background):
    '''
    generate 2x2 contingency table by counting various overlaps

    :param set sample_set: Set of sample genes
    :param set anno_set: Set of annotation genes
    :param BACKGROUND background: Object containing background genes
    :return: 2 dimensional table as described at the top of the program
    '''

    if background is None:
        background = BACKGROUND(anno_set.intersection(sample_set))

    list_anno_overlaps = len(sample_set.intersection(anno_set))
    if list_anno_overlaps == 0:
        return -1
    # alternative code for other specifications
    # background_size = len(sample_set) + len(anno_set) - list_anno_overlaps if background == None else len(background.background_genes)

    background_size = len(background.background_genes)

    list_genome_overlaps = len(sample_set) - list_anno_overlaps
    # alternative code for other specifications
    # genome_anno_overlaps = len(anno_set.intersection(background.background_genes)) if background != None else len(anno_set)
    genome_anno_overlaps = len(anno_set.intersection(background.background_genes))
    genome_only = background_size - genome_anno_overlaps

    return [[list_anno_overlaps, genome_anno_overlaps], [list_genome_overlaps, genome_only]]


def fisher_exact(sample_sets, anno, alpha, background, cpu_count):
    '''
    the fisher exact test

    :param GMT sample_sets: Sample set object
    :param GMT anno: Annotation set object
    :param float alpha: The desired significance level
    :param BACKGROUND background: Object containing background genes
    :param int cpu_count: Specifies the number of cores to be used
    :return: An array of 2 items, the gene rankings and the filtered FDR rankings
    '''

    gene_rankings = []
    # set up multiprocessing inputs
    map_arr = generate_inputs(anno, background)

    for gsid in sample_sets.genesets:

        # each subprocess will append its results to the input_item array
        input_items = multiprocess(gsid, sample_sets, map_arr, fisher_process, cpu_count)
        for input_item in input_items:
            gene_rankings.append(input_item)
    final_rankings = benjamini_hochberg(gene_rankings)
    return [final_rankings, significance_filter(final_rankings, alpha)]

def fisher_process(input_item):
    '''
    fisher sub-method for multiprocessing

    :param list input_item: Object containing input fields
    :return: OverrepResult object
    '''

    cont_table = gen_table(input_item.gene_set, input_item.anno_list, input_item.background)

    # in the case of no overlaps
    if cont_table == -1:
        return OverrepResult(input_item.gsid, len(input_item.gene_set), input_item.anno_id, len(input_item.anno_list),
                             1.0, 0, 0)
    list_anno_overlaps = cont_table[0][0]
    p_value = stats.fisher_exact(cont_table)[1]

    return OverrepResult(input_item.gsid, len(input_item.gene_set), input_item.anno_id, len(input_item.anno_list),
                         p_value, list_anno_overlaps, 0)


def hypergeometric(sample_sets, anno, alpha, background, cpu_count):
    '''
    the hypergeometric test

    :param GMT sample_sets: Sample set object
    :param GMT anno: Annotation set object
    :param float alpha: The desired significance level
    :param BACKGROUND background: Object containing background genes
    :param int cpu_count: Specifies the number of cores to be used
    :return: An array of 2 items, the gene rankings and the filtered FDR rankings
    '''


    gene_rankings = []
    # set up multiprocessing inputs
    map_arr = generate_inputs(anno, background)

    for gsid in sample_sets.genesets:
        # each subprocess will append its results to the input_item array
        input_items = multiprocess(gsid, sample_sets, map_arr, hypergeometric_process, cpu_count)
        for input_item in input_items:
            gene_rankings.append(input_item)
    final_rankings = benjamini_hochberg(gene_rankings)
    return [final_rankings, significance_filter(final_rankings, alpha)]


def hypergeometric_process(input_item):
    '''
    hypergeometric sub-method for multiprocessing

    :param list input_item: Object containing input fields
    :return: OverrepResult object
    '''

    table = gen_table(input_item.gene_set, input_item.anno_list, input_item.background)

    # in the case of no overlaps
    if table == -1:
        return OverrepResult(input_item.gsid, len(input_item.gene_set), input_item.anno_id, len(input_item.anno_list),
                             1.0, 0, 0)
    list_anno_overlaps = table[0][0]
    p_value = stats.hypergeom.sf(table[0][0] - 1, table[0][1] + table[1][1], table[0][1], table[0][0] + table[1][0])

    return OverrepResult(input_item.gsid, len(input_item.gene_set), input_item.anno_id, len(input_item.anno_list),
                         p_value, list_anno_overlaps, 0)


def binomial(sample_sets, anno, alpha, background, cpu_count):
    '''
    the binomial test

    :param GMT sample_sets: Sample set object
    :param GMT anno: Annotation set object
    :param float alpha: The desired significance level
    :param BACKGROUND background: Object containing background genes
    :param int cpu_count: Specifies the number of cores to be used
    :return: An array of 2 items, the gene rankings and the filtered FDR rankings
    '''

    gene_rankings = []
    # set up multiprocessing inputs
    map_arr = generate_inputs(anno, background)
    for gsid in sample_sets.genesets:
        # each subprocess will append its results to the input_item array
        input_items = multiprocess(gsid, sample_sets, map_arr, binomial_process, cpu_count)
        for input_item in input_items:
            gene_rankings.append(input_item)
    final_rankings = benjamini_hochberg(gene_rankings)
    return [final_rankings, significance_filter(final_rankings, alpha)]


def binomial_process(input_item):
    '''
    binomial sub-method for multiprocessing

    :param list input_item: Object containing input fields
    :return: OverrepResult object
    '''

    table = gen_table(input_item.gene_set, input_item.anno_list, input_item.background)

    # in the case of no overlaps
    if table == -1:
        return OverrepResult(input_item.gsid, len(input_item.gene_set), input_item.anno_id, len(input_item.anno_list),
                             1.0, 0, 0)
    list_anno_overlaps = table[0][0]
    p_value = stats.binom_test(list_anno_overlaps, len(input_item.gene_set),
                               float(table[0][1]) / (table[1][1] + table[0][1]))

    return OverrepResult(input_item.gsid, len(input_item.gene_set), input_item.anno_id, len(input_item.anno_list),
                         p_value, list_anno_overlaps, 0)


def chi_squared(sample_sets, anno, alpha, background, cpu_count):
    '''
    the chi squared test

    :param GMT sample_sets: Sample set object
    :param GMT anno: Annotation set object
    :param float alpha: The desired significance level
    :param BACKGROUND background: Object containing background genes
    :param int cpu_count: Specifies the number of cores to be used
    :return: An array of 2 items, the gene rankings and the filtered FDR rankings
    '''

    gene_rankings = []
    # set up multiprocessing inputs
    map_arr = generate_inputs(anno, background)
    for gsid in sample_sets.genesets:
        # each subprocess will append its results to the input_item array
        input_items = multiprocess(gsid, sample_sets, map_arr, chi_process, cpu_count)
        for input_item in input_items:
            gene_rankings.append(input_item)

    final_rankings = benjamini_hochberg(gene_rankings)
    return [final_rankings, significance_filter(final_rankings, alpha)]

def chi_process(input_item):
    '''
    chi squared sub-method for multiprocessing

    :param list input_item: Object containing input fields
    :return: OverrepResult object
    '''

    table = gen_table(input_item.gene_set, input_item.anno_list, input_item.background)

    # in the case of no overlaps
    if table == -1:
        return OverrepResult(input_item.gsid, len(input_item.gene_set), input_item.anno_id, len(input_item.anno_list),
                             1.0, 0, 0)
    list_anno_overlaps = table[0][0]
    p_value = stats.chisquare(table)[1][0]
    return OverrepResult(input_item.gsid, len(input_item.gene_set), input_item.anno_id, len(input_item.anno_list),
                         p_value, list_anno_overlaps, 0)


def benjamini_hochberg(gene_rankings):
    '''
    benjamini hochberg FDR correction for multiple hypothesis testing

    :param list gene_rankings: Object containing OverrepResult objects
    :return: A p value adjusted array of OverrepResult objects
    '''

    output = []
    # sort rankings before correcting p values
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line.p_value))

    prev_bh_value = 0
    for i, OR_Result in enumerate(gene_rankings):

        # standard BH correction based on ranking
        bh_value = OR_Result.p_value * len(gene_rankings) / float(i + 1)
        bh_value = min(bh_value, 1)

        # to preserve monotonicity
        # ensures that less significant items do not have a higher FDR than more significant items
        if bh_value < prev_bh_value:
            output[i - 1].FDR = bh_value

        OR_Result.FDR = bh_value
        output.append(OR_Result)
        prev_bh_value = bh_value

    return output


def significance_filter(gene_rankings, alpha):
    '''
    filters out significant items under the alpha level

    :param list gene_rankings: Object containing OverrepResult objects
    :param float alpha: The desired significance level
    :return: A filtered array of OverrepResult objects
    '''

    significant_values = []
    for i, OR_Result in enumerate(gene_rankings):
        if OR_Result.FDR <= alpha:
            significant_values.append(OR_Result)

    return significant_values


if __name__ == '__main__':
    from argparse import ArgumentParser

    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, version="%prog dev-unreleased")

    parser.add_argument(
        "-n",
        "--name of test to run",
        dest="test_name",
        help="options --> (fisher_exact, hypergeometric, binomial, chi_squared)",
        metavar="STRING",
        required=True
    )

    parser.add_argument(
        "-g",
        "--gene-sets file",
        dest="gene_sets",
        help="gene sets file for comparision",
        metavar=".gmt FILE",
        required=True
    )

    parser.add_argument(
        "-a",
        "--annotations-file",
        dest="annotation_list",
        help="annotation file",
        metavar=".gmt FILE",
        required=True,
        default=None
    )
    parser.add_argument(
        "-b",
        "--background-gene file",
        dest="background_list",
        help="background gene list file for comparision (DEFAULT None)",
        metavar=".txt FILE"
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

    # perform a fisher_exact test with a console printout and including all values
    test = OverrepTest()
    test.run()
