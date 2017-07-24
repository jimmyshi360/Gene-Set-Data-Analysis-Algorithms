import multiprocessing
import time

from flib.core.gmt import GMT
from scipy import stats

from overrep_ouput_writer import OUT
from utilities.background import BACKGROUND

'''
FOR PARSER OTIONS:
-g = GENE LIST (.gmt)
-a = ANNOTATIONS LIST (.gmt)
-b = BACKGROUND LIST (.gmt)
-c = CLUSTER NUMBER (integer)
-l = CLUSTER LIST (.mat)
-o = OUTPUT (.txt)
-r = FALSE DISCOVERY RATE (float)

2x2 Contingency Table for Chi Squared and Fisher Exact Tests

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

class EnrichmentResult:
    def __init__(self, gsid, gs1_ngenes, go_id, gs2_ngenes,p_value, overlaps, FDR):
        self.gsid=gsid
        self.gs1_ngenes=gs1_ngenes
        self.go_id = go_id
        self.gs2_ngenes=gs2_ngenes
        self.p_value = p_value
        self.overlaps=overlaps
        self.FDR=FDR

    def set_FDR(self,new_FDR):
        self.FDR=new_FDR

class InputItem:
    def __init__(self, go_id, gene_list, background,gsid, gene_set):
        self.go_id = go_id
        self.gene_list=gene_list
        self.background=background
        self.gsid=gsid
        self.gene_set=gene_set


p = None
args=None

# generates a list of inputs to be mapped to several processors
# each parameter will be an array of terms
def generate_inputs(anno, background):

    input_arr = []
    for go_id in anno.genesets:
        input_arr.append([go_id, anno.genesets[go_id], background])
    return input_arr


# completes the input array for each sepcific sample gene set and executes a multiprocess mapping inputs to all GO gene sets
def multiprocess(gsid, sample, map_arr, method):
    input_arr=[]
    for i, row in enumerate(map_arr):
        go_id=row[0]
        gene_list=row[1]
        background=row[2]
        next_item = InputItem(go_id, gene_list, background, gsid, sample.genesets[gsid])
        input_arr.append(next_item)

    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = p.map(method, input_arr)
    p.close()
    p.join()

    return results


# generate contingency table
def gen_table(sample_set, anno_set, background):

    if background is None:
        background=BACKGROUND(anno_set.intersection(sample_set))

    list_anno_overlaps = len(sample_set.intersection(anno_set))
    if list_anno_overlaps==0:
        return 0
    #background_size = len(sample_set) + len(anno_set) - list_anno_overlaps if background == None else len(background.background_genes)
    background_size =len(background.background_genes)

    list_genome_overlaps = len(sample_set) - list_anno_overlaps
    #genome_anno_overlaps = len(anno_set.intersection(background.background_genes)) if background != None else len(anno_set)
    genome_anno_overlaps = len(anno_set.intersection(background.background_genes))
    genome_only = background_size - genome_anno_overlaps

    return [[list_anno_overlaps, genome_anno_overlaps], [list_genome_overlaps, genome_only]]


# fisher exact test on 2x2 contigency tables
def fisher_exact(sample, anno, alpha, background):

    gene_rankings = []
    # contruct and analyze contingency tables
    map_arr = generate_inputs(anno, background)

    for gsid in sample.genesets:
        items = multiprocess(gsid, sample, map_arr, fisher_process)
        for item in items:
            gene_rankings.append(item)
    final_rankings = benjamini_hochberg(gene_rankings)
    return [final_rankings, significance_filter(final_rankings, alpha)]


# fisher sub-method for multiprocessing
def fisher_process(input_item):

    table = gen_table(input_item.gene_set,input_item.gene_list, input_item.background)

    # in the case of no overlaps
    if table==0:
        return EnrichmentResult(input_item.gsid, len(input_item.gene_set), input_item.go_id, len(input_item.gene_list), 1.0, 0,0)
    list_anno_overlaps = table[0][0]
    p_value = stats.fisher_exact(table)[1]

    return EnrichmentResult(input_item.gsid, len(input_item.gene_set), input_item.go_id, len(input_item.gene_list), p_value, list_anno_overlaps,0)


# hypergeometric test on 2x2 contigency tables
def hypergeometric(sample, anno, alpha,background):
    gene_rankings = []
    map_arr = generate_inputs(anno, background)

    for gsid in sample.genesets:
        items = multiprocess(gsid, sample, map_arr, hypergeometric_process)
        for item in items:
            gene_rankings.append(item)
    final_rankings = benjamini_hochberg(gene_rankings)
    return [final_rankings, significance_filter(final_rankings, alpha)]

# hypergeometric sub-method for multiprocessing, m_arr contains the parameters
def hypergeometric_process(input_item):
    table = gen_table(input_item.gene_set, input_item.gene_list, input_item.background)

    # in the case of no overlaps
    if table == 0:
        return EnrichmentResult(input_item.gsid, len(input_item.gene_set), input_item.go_id, len(input_item.gene_list),
                                1.0, 0, 0)
    list_anno_overlaps = table[0][0]
    p_value = stats.hypergeom.sf(table[0][0] - 1, table[0][1] + table[1][1], table[0][1], table[0][0] + table[1][0])

    return EnrichmentResult(input_item.gsid, len(input_item.gene_set), input_item.go_id, len(input_item.gene_list), p_value, list_anno_overlaps,0)


# binomial test on 2x2 contingency tables
def binomial(sample, anno, alpha,background):
    gene_rankings = []
    map_arr = generate_inputs(anno, background)
    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        items = multiprocess(gsid, sample, map_arr, binomial_process)
        for item in items:
            gene_rankings.append(item)
    final_rankings = benjamini_hochberg(gene_rankings)
    return [final_rankings, significance_filter(final_rankings, alpha)]


# binomial sub-method for multiprocessing
def binomial_process(input_item):
    table = gen_table(input_item.gene_set, input_item.gene_list, input_item.background)
    list_anno_overlaps = table[0][0]

    # in the case of no overlaps
    if table == 0:
        return EnrichmentResult(input_item.gsid, len(input_item.gene_set), input_item.go_id, len(input_item.gene_list),
                                1.0, 0, 0)
    list_anno_overlaps = table[0][0]
    p_value = stats.binom_test(list_anno_overlaps, len(input_item.gene_set), float(table[0][1]) / (table[1][1] + table[0][1]))

    return EnrichmentResult(input_item.gsid, len(input_item.gene_set), input_item.go_id, len(input_item.gene_list),p_value, list_anno_overlaps,0)


# chi squared test on 2x2 contigency tables
def chi_squared(sample, anno, alpha,background):
    gene_rankings = []
    map_arr = generate_inputs(anno, background)
    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        items = multiprocess(gsid, sample, map_arr, chi_process)
        for item in items:
            gene_rankings.append(item)

    final_rankings=benjamini_hochberg(gene_rankings)
    return [final_rankings, significance_filter(final_rankings,alpha)]

# chi squared sub-method for multiprocessing
def chi_process(input_item):
    table = gen_table(input_item.gene_set, input_item.gene_list, input_item.background)


    # in the case of no overlaps
    if  table == 0:
        return EnrichmentResult(input_item.gsid, len(input_item.gene_set), input_item.go_id, len(input_item.gene_list),
                                1.0, 0, 0)
    list_anno_overlaps = table[0][0]
    p_value = stats.chisquare(table)[1][0]
    return EnrichmentResult(input_item.gsid, len(input_item.gene_set), input_item.go_id, len(input_item.gene_list),p_value, list_anno_overlaps,0)

# wrapper method for integrating parsers and all overrep tests
def over_rep_test(test_name, print_option,background=None):

    sample = GMT(args.gene_list)
    anno = GMT(args.annotation_list)
    if args.background_list is not None:
        background = BACKGROUND([],args.background_list)

    if test_name == "fisher_exact":
        rankings = fisher_exact(sample, anno, args.rate,background)
    elif test_name == "chi_squared":
        rankings = chi_squared(sample, anno, args.rate, background)
    elif test_name == "binomial":
        rankings = binomial(sample, anno,args.rate, background)
    elif test_name == "hypergeometric":
        rankings = hypergeometric(sample, anno, args.rate,background)

    # prints out the rankings and significant values

    return OUT(rankings[0],rankings[1], args.output).printout(print_option,False)

# FDR correction for multiple hypothesis testing
def benjamini_hochberg(gene_rankings):

    output=[]
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line.p_value))

    prev_bh_value=0
    for i, E_Result in enumerate(gene_rankings):

        bh_value = E_Result.p_value * len(gene_rankings) / float(i + 1)
        bh_value = min(bh_value, 1)

        #to preserve monotonicity
        if bh_value<prev_bh_value:
            output[i-1].set_FDR(bh_value)

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
        "-g",
        "--gene-list file",
        dest="gene_list",
        help="gene list file for comparision",
        metavar=".gmt FILE",
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
        "-b",
        "--background-gene file",
        dest="background_list",
        help="background gene list file for comparision (OPTIONAL)",
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

    args = parser.parse_args()

    over_rep_test("fisher_exact", True)
