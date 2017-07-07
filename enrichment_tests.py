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

______________________________Gene List_______________All Observed Genes___________
In Annotations     |    list_anno_overlaps    |       genome_anno_overlaps         |
------------------------------------------------------------------------------------
Not in Annotations |    list_genome_overlaps  |          genome_only               |
------------------------------------------------------------------------------------
'''

p=None

# fisher exact test on 2x2 contigency tables
def fisher_exact(sample, anno, background):
    background_size = 0 if background == None else len(background.background_genes)
    gene_rankings=[]
    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        for go_id in anno.genesets:
            current_set=sample.genesets[gsid]
            list_anno_overlaps = len(current_set.intersection(anno.genesets[go_id]))
            list_genome_overlaps = len(current_set) - list_anno_overlaps
            genome_anno_overlaps = len(anno.genesets[go_id]) - list_anno_overlaps
            genome_only = background_size - list_anno_overlaps - list_genome_overlaps - genome_anno_overlaps if background != None else 0
            p_value = stats.fisher_exact([[list_anno_overlaps, genome_anno_overlaps], [list_genome_overlaps, genome_only]])[1]

            gene_rankings.append([gsid, go_id, p_value,  list_anno_overlaps])

    return gene_rankings

# hypergeometric test on 2x2 contigency tables
def hypergeometric(sample, anno, background):
    background_size = 0 if background == None else len(background.background_genes)
    gene_rankings = []

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        for go_id in anno.genesets:
            list_anno_overlaps = len(sample.genesets[gsid].intersection(anno.genesets[go_id]))
            list_genome_overlaps = len(sample.genesets[gsid]) - list_anno_overlaps
            genome_anno_overlaps = len(anno.genesets[go_id]) - list_anno_overlaps
            genome_only = background_size - list_anno_overlaps - list_genome_overlaps - genome_anno_overlaps if background != None else 0
            p_value =stats.hypergeom.sf(list_anno_overlaps-1, genome_anno_overlaps+genome_only,  genome_anno_overlaps, list_anno_overlaps+list_genome_overlaps)
            gene_rankings.append([gsid, go_id, p_value,  list_anno_overlaps])
    return gene_rankings


#binomial test on 2x2 contingency tables
def binomial(sample, anno, background):
    background_size = 0 if background == None else len(background.background_genes)
    gene_rankings = []

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        for go_id in anno.genesets:
            list_anno_overlaps = len(sample.genesets[gsid].intersection(anno.genesets[go_id]))
            list_genome_overlaps = len(sample.genesets[gsid]) - list_anno_overlaps
            genome_anno_overlaps = len(anno.genesets[go_id]) - list_anno_overlaps
            genome_only = background_size if background != None else genome_anno_overlaps
            p_value = stats.binom_test(list_anno_overlaps, len(sample.genesets[gsid]), float(genome_anno_overlaps)/genome_only)
            gene_rankings.append([gsid, go_id, p_value,  list_anno_overlaps])

    return gene_rankings

# chi squared test on 2x2 contigency tables
def chi_squared(sample, anno, background):
    background_size = 0 if background == None else len(background.background_genes)
    gene_rankings = []

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        for go_id in anno.genesets:
            list_anno_overlaps = len(sample.genesets[gsid].intersection(anno.genesets[go_id]))
            list_genome_overlaps = len(sample.genesets[gsid]) - list_anno_overlaps
            genome_anno_overlaps = len(anno.genesets[go_id]) - list_anno_overlaps
            genome_only = background_size - list_anno_overlaps - list_genome_overlaps - genome_anno_overlaps if background != None else 0
            p_value = stats.chisquare([[list_anno_overlaps, genome_anno_overlaps], [list_genome_overlaps, genome_only]])[1][0]
            gene_rankings.append([gsid, go_id, p_value,  list_anno_overlaps])

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
        background = BACKGROUND(p.args.background_list)
        use_parsers = True
    else:
        sample = GMT(sample)
        anno = GMT(anno)
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
    def test_binomial(self):
        self.assertEqual(len(over_rep_test("binomial", False, "GMT.gmt", "GO.gmt", "BACKGROUND.txt", 0.05, "OUTPUT.txt")),
                         1264)
        self.assertEqual(over_rep_test("binomial", False,"GMT.gmt", "GO.gmt", "BACKGROUND.txt", 0.05, "OUTPUT.txt")[3],
                         ["0", "17", "GO:0090083", "6", "0", "1.0", "0.0297153710663"])

    def test_fisher(self):
        self.assertEqual(len(over_rep_test("fisher_exact", False,"GMT.gmt", "GO.gmt", "BACKGROUND.txt", 0.05, "OUTPUT.txt")),
                         1263)
        self.assertEqual(over_rep_test("fisher_exact", False,"GMT.gmt", "GO.gmt", "BACKGROUND.txt", 0.05, "OUTPUT.txt")[3],
                         ["0", "17", "GO:0090083", "6", "0", "1.0", "0.0303094331668"])

    def test_chi_squared(self):
        self.assertEqual(
            len(over_rep_test("chi_squared", False,"GMT.gmt", "GO.gmt", "BACKGROUND.txt", 0.05, "OUTPUT.txt")), 40106)
        self.assertEqual(over_rep_test("chi_squared", False,"GMT.gmt", "GO.gmt", "BACKGROUND.txt", 0.05, "OUTPUT.txt")[3],
                         ["0", "17", "GO:0090344", "6", "0", "3.73798184017e-05", "0.045758367941"])

    def test_hypergeometric(self):
        self.assertEqual(
            len(over_rep_test("hypergeometric", False,"GMT.gmt", "GO.gmt", "BACKGROUND.txt", 0.05, "OUTPUT.txt")), 1268)
        self.assertEqual(over_rep_test("hypergeometric", False,"GMT.gmt", "GO.gmt", "BACKGROUND.txt", 0.05, "OUTPUT.txt")[3],
                         ["0", "17", "GO:0090083", "6", "0", "1.0", "0.0295011396605"])

if __name__ == '__main__':

    #choose fisher_exact, chi_squared, hypergeometric, or binomial
    print("TESTING")
suite = unittest.TestLoader().loadTestsFromTestCase(TestStattests)
unittest.TextTestRunner(verbosity=2).run(suite)

