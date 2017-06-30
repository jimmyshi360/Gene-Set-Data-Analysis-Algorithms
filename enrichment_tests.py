from scipy import stats
from flib.core.gmt import GMT
from output_generator import OUT
from background import BACKGROUND
from mat import MAT
import math
import scipy.stats
import numpy as np
from parser_adder import Parsers

'''
FOR PARSER CHOICES:
-g = GENE LIST (.gmt)
-a = ANNOTATIONS LIST (.gmt)
-b = BACKGROUND LIST (.gmt)
-c = CLUSTER NUMBER (integer)
-l = CLUSTER LIST (.mat)
-o = OUTPUT (.txt)
-r = FALSE DISCOVERY RATE (float)
'''

'''2x2 Contingency Table for Chi Squared and Fisher Exact Tests
______________________________Gene List_______________All Observed Genes___________
In Annotations     |                          |                                    |
------------------------------------------------------------------------------------
Not in Annotations |  list_genome_overlaps    |                                    |
------------------------------------------------------------------------------------
'''

p=None

# fisher exact test on 2x2 contigency tables
def fisher_exact():

    p = Parsers("-g -a -b -o -r")
    sample = GMT(p.args.gene_list)
    anno = GMT(p.args.annotation_list)
    background = BACKGROUND(p.args.background_list)

    gene_rankings = []
    anno_genes = anno._genes
    background_size = 0 if background == None else len(background.background_genes)

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        list_anno_overlaps = len(sample.genesets[gsid].intersection(anno_genes))
        list_genome_overlaps = len(sample.genesets[gsid]) - list_anno_overlaps
        genome_anno_overlaps = len(anno_genes) - list_anno_overlaps
        genome_only = background_size - list_anno_overlaps - list_genome_overlaps - genome_anno_overlaps if background != None else 0
        p_value = stats.fisher_exact([[list_anno_overlaps, genome_anno_overlaps], [list_genome_overlaps, genome_only]])[1]
        gene_rankings.append([p_value, gsid])

    # prints out the rankings and significant values
    OUT(gene_rankings, p.args.output, p.args.rate).printout()

# chi squared test on 2x2 contigency tables
def chi_square():

    p = Parsers("-g -a -b -o -r")
    sample = GMT(p.args.gene_list)
    anno = GMT(p.args.annotation_list)
    background = BACKGROUND(p.args.background_list)

    gene_rankings = []
    anno_genes = anno._genes
    background_size = 0 if background == None else len(background.background_genes)

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        list_anno_overlaps = len(sample.genesets[gsid].intersection(anno_genes))
        list_genome_overlaps = len(sample.genesets[gsid]) - list_anno_overlaps
        genome_anno_overlaps = len(anno_genes) - list_anno_overlaps
        genome_only = background_size - list_anno_overlaps - list_genome_overlaps - genome_anno_overlaps if background != None else 0
        p_value = stats.chisquare([[list_anno_overlaps, genome_anno_overlaps], [list_genome_overlaps, genome_only]],1)[1][0]
        gene_rankings.append([p_value, gsid])

    # prints out the rankings and significant values
    OUT(gene_rankings, p.args.output, p.args.rate).printout()

#wilcoxon rank sum test, compares an input list of genesets versus scores between two experimental groups
def wilcoxon():
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
        p_value = scipy.stats.ranksums(score_arr, total_score_list)

        gene_rankings.append([p_value[1], gsid])
        gene_rankings = sorted(gene_rankings, key=lambda line: float(line[0]))

    # prints out the rankings and significant values
    OUT(gene_rankings, p.args.output, p.args.rate).printout()

#parametric analysis gene enrichment test, compares an input list of genesets versus scores between two experimental groups
def page():

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
    OUT(gene_rankings, p.args.output, p.args.rate).printout()

#binomial test on 2x2 contingency tables
def binomial():

    p = Parsers("-g -a -b -o -r")
    sample = GMT(p.args.gene_list)
    anno = GMT(p.args.annotation_list)
    background = BACKGROUND(p.args.background_list)

    gene_rankings = []
    anno_genes = anno._genes
    background_size = 0 if background == None else len(background.background_genes)

    # contruct and analyze contingency tables
    for gsid in sample.genesets:
        list_anno_overlaps = len(sample.genesets[gsid].intersection(anno_genes))
        list_genome_overlaps = len(sample.genesets[gsid]) - list_anno_overlaps
        genome_anno_overlaps = len(anno_genes) - list_anno_overlaps
        genome_only = background_size if background != None else genome_anno_overlaps
        p_value = stats.binom_test(list_anno_overlaps, len(sample.genesets[gsid]), float(genome_anno_overlaps)/genome_only)
        gene_rankings.append([p_value, gsid])

    # prints out the rankings and significant values
    OUT(gene_rankings, p.args.output, p.args.rate).printout()

#switch case for choosoing the test type
def test_switch(x):
    return {
        1 : fisher_exact(),
        2: page(),
        3: binomial(),
        4: chi_square(),
        5: wilcoxon()
    }[x]

if __name__ == '__main__':

    #print("Please choose which test to run\n1) Fisher Exact\n2) PAGE\n3) Binomial\n4) chi_square()\n5) wilcoxon()")
    #test_switch(int(input()))

    fisher_exact(),
    #page(),
    #binomial(),
    #chi_square(),
    #wilcoxon()