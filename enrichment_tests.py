from scipy import stats
from flib.core.gmt import GMT
from anno_ouput_writer import OUT
from background import BACKGROUND
from parsers import Parsers
from mat import MAT
import math
import sys
import scipy.stats
import numpy as np
from collections import defaultdict
from itertools import islice
import matplotlib.pyplot as plt
import unittest

import time
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
'''

p=None

#generate contingency table
def gen_table(sample_set,anno_set,background):
    list_anno_overlaps = len(sample_set.intersection(anno_set))
    background_size = len(sample_set)+len(anno_set)-list_anno_overlaps if background == None else len(background.background_genes)
    list_genome_overlaps = len(sample_set) - list_anno_overlaps
    genome_anno_overlaps = len(anno_set.intersection(background.background_genes)) if background!= None else len(anno_set)
    genome_only = background_size - genome_anno_overlaps
    return [[list_anno_overlaps, genome_anno_overlaps],[list_genome_overlaps, genome_only]]

# fisher exact test on 2x2 contigency tables
def fisher_exact(sample, anno, background):

    gene_rankings=[]
    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        for go_id in anno.genesets:
            table=gen_table(sample.genesets[gsid],anno.genesets[go_id],background)
            p_value = stats.fisher_exact(table)[1]
            gene_rankings.append([gsid, go_id, p_value,  table[0][0]])

    return gene_rankings

# hypergeometric test on 2x2 contigency tables
def hypergeometric(sample, anno, background):

    gene_rankings = []

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        for go_id in anno.genesets:
            table = gen_table(sample.genesets[gsid], anno.genesets[go_id], background)
            p_value =stats.hypergeom.sf(table[0][0]-1,table[0][1]+table[1][1],  table[0][1], table[0][0]+table[1][0])
            gene_rankings.append([gsid, go_id, p_value,  table[0][0]])
    return gene_rankings


#binomial test on 2x2 contingency tables
def binomial(sample, anno, background):

    gene_rankings = []

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        for go_id in anno.genesets:
            table = gen_table(sample.genesets[gsid], anno.genesets[go_id], background)
            p_value = stats.binom_test(table[0][0], len(sample.genesets[gsid]), float(table[0][1])/(table[1][1]+table[0][1]))
            gene_rankings.append([gsid, go_id, p_value,  table[0][0]])

    return gene_rankings

# chi squared test on 2x2 contigency tables
def chi_squared(sample, anno, background):

    gene_rankings = []

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        for go_id in anno.genesets:
            table = gen_table(sample.genesets[gsid], anno.genesets[go_id], background)
            p_value = stats.chisquare(table)[1][0]
            gene_rankings.append([gsid, go_id, p_value,  table[0][0]])

    return gene_rankings

#wilcoxon rank sum test, compares an input list of genesets versus scores between two experimental groups
def wilcoxon():
    print("\nWILCOXON")
    p = Parsers("-g -c -o -l -r")
    gmt = GMT(p.args.gene_list)
    mat = MAT(p.args.cluster_list)
    cluster=p.args.cluster_number

    score_arr = []
    gene_rankings = []

    for gsid in gmt.genesets:

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
    return OUT(gene_rankings, p.args.output, p.args.rate).printout()

#parametric analysis gene enrichment test, compares an input list of genesets versus scores between two experimental groups
def page():
    print("\nPAGE")
    p=Parsers("-g -c -o -l -r")
    gmt = GMT(p.args.gene_list)
    mat = MAT(p.args.cluster_list)
    cluster=p.args.cluster_number
    gene_rankings = []
    gene_mean = 0
    gene_sd = 0
    geneset_mean = 0
    geneset_size = 0

    score_arr = []
    #calculate value related to the entire cluster
    for i in mat.matrix.keys():
        score_arr.append(list(mat.matrix[i])[cluster])
    score_arr = np.array(score_arr).astype(np.float)
    gene_mean=np.mean(score_arr)
    gene_sd = np.std(score_arr)

    #calculate p values based on mean, standard deviation and list sizes
    for gsid in gmt.genesets:
        geneset_size = len(gmt.genesets[gsid])
        score_arr = []
        #for each gene set, calculate values
        for gene in gmt.genesets[gsid]:
             row_arr = list(mat.matrix[gene])
             score_arr.append(row_arr[cluster])
        score_arr = np.array(score_arr).astype(np.float)
        gene_set_mean = np.mean(score_arr)
        z_score = (gene_mean-geneset_mean) * math.sqrt(geneset_size) / gene_sd
        p_value = stats.norm.sf(abs(z_score))
        gene_rankings.append([p_value, gsid])

    #prints out the rankings and significant values
        return OUT(gene_rankings, p.args.output, p.args.rate).printout()

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
        if background!=None:
            background = BACKGROUND(background)

    if test_name=="fisher_exact":
        gene_rankings=fisher_exact(sample,anno,background)
    elif test_name=="chi_squared":
        gene_rankings=chi_squared(sample,anno,background)
    elif test_name=="binomial":
        gene_rankings=binomial(sample,anno,background)
    elif test_name=="hypergeometric":
        gene_rankings=hypergeometric(sample,anno,background)

    # prints out the rankings and significant values
    if use_parsers:
        return OUT(gene_rankings, p.args.output, p.args.rate, sample, anno).printout(print_option)
    else:
        return OUT(gene_rankings, output, rate, sample, anno).printout(print_option)


class TestStattests(unittest.TestCase):


    def test_gen_table(self):
        anno = GMT("unittest_files\go.gmt").genesets['0']
        sample = GMT("unittest_files\gmt.gmt").genesets['0']
        background = BACKGROUND("unittest_files\\background.txt")

        self.assertEqual(gen_table(sample, anno, background)[0][0], 3)
        self.assertEqual(gen_table(sample, anno, background)[0][1], 40)
        self.assertEqual(gen_table(sample, anno, background)[1][0], 297)
        self.assertEqual(gen_table(sample, anno, background)[1][1], 19960)

    def test_binomial(self):
        anno = "unittest_files\go.gmt"
        sample = "unittest_files\gmt.gmt"
        background = "unittest_files\\background.txt"
        output= "unittest_files\output.txt"
        self.assertEqual(over_rep_test("binomial", False, sample, anno, background, 0.05, output)[0][5], "0.0229768421702")

    def test_fisher(self):

        anno = "unittest_files\go.gmt"
        sample =  "unittest_files\gmt.gmt"
        background = "unittest_files\\background.txt"
        output = "unittest_files\output.txt"
        self.assertEqual(over_rep_test("fisher_exact", False,sample, anno, background, 0.05, output)[0][5],"0.0255246814673")


    def test_chi_squared(self):

        anno = "unittest_files\go.gmt"
        sample = "unittest_files\gmt.gmt"
        background = "unittest_files\\background.txt"
        output = "unittest_files\output.txt"
        self.assertEqual(over_rep_test("chi_squared", False,sample, anno, background, 0.05, output)[0][5], "1.27701446634e-64")


    def test_hypergeometric(self):

        anno = "unittest_files\go.gmt"
        sample = "unittest_files\gmt.gmt"
        background = "unittest_files\\background.txt"
        output = "unittest_files\output.txt"

        self.assertEqual(
            over_rep_test("hypergeometric", False,sample, anno, background, 0.05, output)[0][5], "0.0219349067622")



if __name__ == '__main__':

    #choose fisher_exact, chi_squared, hypergeometric, or binomial
    anno = "unittest_files\go.gmt"
    sample = "unittest_files\gmt.gmt"
    background = "unittest_files\\background.txt"
    output = "unittest_files\output.txt"
    over_rep_test("binomial", True, sample, anno, background, 0.05, output)

suite = unittest.TestLoader().loadTestsFromTestCase(TestStattests)
unittest.TextTestRunner(verbosity=2).run(suite)

