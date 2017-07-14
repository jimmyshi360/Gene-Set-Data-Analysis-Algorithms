import multiprocessing

from flib.core.gmt import GMT
from scipy import stats

from overrep_ouput_writer import OUT
from utilities.background import BACKGROUND
from utilities.parsers import Parsers

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

p = None


# generates a list of inputs to be mapped to several processors
# each parameter will be an array of terms
def build_inputs(anno, background):
    input_arr = []
    for go_id in anno.genesets:
        input_arr.append([go_id, anno.genesets[go_id], background])
    return input_arr


# completes the input array for each sepcific sample gene set and executes a multiprocess mapping inputs to all GO gene sets
def multiprocess(gsid, sample, map_arr, method):
    map_arr_copy = list(map_arr)
    for i in range(0, len(map_arr)):
        map_arr_copy[i] = [map_arr_copy[i][0], map_arr_copy[i][1], map_arr_copy[i][2], gsid, sample.genesets[gsid]]

    p = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = p.map(method, map_arr_copy)
    p.close()
    p.join()

    return results


# generate contingency table
def gen_table(sample_set, anno_set, background):
    list_anno_overlaps = len(sample_set.intersection(anno_set))
    background_size = len(sample_set) + len(anno_set) - list_anno_overlaps if background == None else len(
        background.background_genes)
    list_genome_overlaps = len(sample_set) - list_anno_overlaps
    genome_anno_overlaps = len(anno_set.intersection(background.background_genes)) if background != None else len(
        anno_set)
    genome_only = background_size - genome_anno_overlaps
    return [[list_anno_overlaps, genome_anno_overlaps], [list_genome_overlaps, genome_only]]


# fisher exact test on 2x2 contigency tables
def fisher_exact(sample, anno, background):
    gene_rankings = []
    # contruct and analyze contingency tables
    map_arr = build_inputs(anno, background)

    for gsid in sample.genesets:
        items = multiprocess(gsid, sample, map_arr, fisher_process)
        for i in range(0, len(items)):
            gene_rankings.append(items[i])

    return gene_rankings


# fisher sub-method for multiprocessing
def fisher_process(m_arr):
    table = gen_table(m_arr[4], m_arr[1], m_arr[2])
    p_value = stats.fisher_exact(table)[1]
    return [m_arr[3], m_arr[0], p_value, table[0][0]]


# hypergeometric test on 2x2 contigency tables
def hypergeometric(sample, anno, background):
    gene_rankings = []
    map_arr = build_inputs(anno, background)

    for gsid in sample.genesets:
        items = multiprocess(gsid, sample, map_arr, hypergeometric_process)
        for i in range(0, len(items)):
            gene_rankings.append(items[i])
    return gene_rankings

# hypergeometric sub-method for multiprocessing
def hypergeometric_process(m_arr):
    table = gen_table(m_arr[4], m_arr[1], m_arr[2])
    p_value = stats.hypergeom.sf(table[0][0] - 1, table[0][1] + table[1][1], table[0][1], table[0][0] + table[1][0])
    return [m_arr[3], m_arr[0], p_value, table[0][0]]


# binomial test on 2x2 contingency tables
def binomial(sample, anno, background):
    gene_rankings = []
    map_arr = build_inputs(anno, background)
    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        items = multiprocess(gsid, sample, map_arr, binomial_process)
        for i in range(0, len(items)):
            gene_rankings.append(items[i])
    return gene_rankings

# binomial sub-method for multiprocessing
def binomial_process(m_arr):
    table = gen_table(m_arr[4], m_arr[1], m_arr[2])
    p_value = stats.binom_test(table[0][0], len(m_arr[4]), float(table[0][1]) / (table[1][1] + table[0][1]))
    return [m_arr[3], m_arr[0], p_value, table[0][0]]


# chi squared test on 2x2 contigency tables
def chi_squared(sample, anno, background):
    gene_rankings = []
    map_arr = build_inputs(anno, background)
    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        items = multiprocess(gsid, sample, map_arr, chi_process)
        for i in range(0, len(items)):
            gene_rankings.append(items[i])
    return gene_rankings

# chi squared sub-method for multiprocessing
def chi_process(m_arr):
    table = gen_table(m_arr[4], m_arr[1], m_arr[2])
    p_value = stats.chisquare(table)[1][0]
    return [m_arr[3], m_arr[0], p_value, table[0][0]]

# wrapper method for integrating parsers and all overrep tests
def over_rep_test(test_name, print_option, sample=None, anno=None, background=None, rate=None, output=None):
    use_parsers = False
    if sample == None:
        p = Parsers("-g -a -b -o -r")
        sample = GMT(p.args.gene_list)
        anno = GMT(p.args.annotation_list)
        if p.args.background_list != None:
            background = BACKGROUND(p.args.background_list)
        use_parsers = True
    else:
        sample = GMT(sample)
        anno = GMT(anno)
        if background != None:
            background = BACKGROUND(background)

    if test_name == "fisher_exact":
        gene_rankings = fisher_exact(sample, anno, background)
    elif test_name == "chi_squared":
        gene_rankings = chi_squared(sample, anno, background)
    elif test_name == "binomial":
        gene_rankings = binomial(sample, anno, background)
    elif test_name == "hypergeometric":
        gene_rankings = hypergeometric(sample, anno, background)

    # prints out the rankings and significant values
    if use_parsers:
        return OUT(gene_rankings, p.args.output, p.args.rate, sample, anno).printout(print_option)
    else:
        return OUT(gene_rankings, output, rate, sample, anno).printout(print_option)


if __name__ == '__main__':
    # choose fisher_exact, chi_squared, hypergeometric, or binomial
    anno = "test_files\GO.gmt"
    sample = "test_files\GMT.gmt"
    background = "test_files\BACKGROUND.txt"
    output = "test_files\OUTPUT.txt"
    over_rep_test("fisher_exact", True, sample, anno, background, 0.05, output)
