import math
import sys

import numpy as np
from flib.core.gmt import GMT
from scipy import stats
from work_in_progress.enrichment_output_writer import E_OUT

from overrep_ouput_writer import OUT
from utilities.mat import MAT
from utilities.parsers import Parsers


def enrichment_score(anno, cluster, mat, weight):

    N = len(anno)

    set_score_sum=0
    for id in mat.keys():
        if len(list(mat[id]))!=0:
            set_score_sum+=(float(list(mat[id])[cluster]))**weight

    max_ES = -sys.maxint
    min_ES = sys.maxint
    net_sum = 0
    sum_arr = []
    for i in mat.keys():
        if i in anno:
            if len(list(mat[id])) != 0:
                net_sum += (float(list(mat[i])[cluster]))**weight/set_score_sum
        else:
            net_sum += -1.0/(len(mat)-len(anno))
        if net_sum > max_ES:
            max_ES = net_sum
        if net_sum < min_ES:
            min_ES = net_sum
        sum_arr.append(net_sum)

    return max_ES

def permute(mat, anno):
    permuted_arr=[]
    for i in range(0,50):
        anno_permute = np.random.permutation(list(anno))
        permuted_arr.append(anno_permute[0:len(mat)])
    return permuted_arr

def es_distr(mat,cluster, permuted_arr):
    es_arr=[]
    for i in range(0,len(permuted_arr)):
        es_arr.append(enrichment_score(permuted_arr[i],cluster,mat,1))
    return es_arr

def gsea(mat,anno, cluster):

    gene_rankings=[]

    for go_id in anno.genesets:
        es=enrichment_score(anno.genesets[go_id], cluster, mat.matrix, 1)
        permuted_arr = permute(mat.matrix, anno.genesets[go_id])
        es_arr=es_distr(mat.matrix, cluster, permuted_arr)
        nes=normalize_score(es,es_arr)
        nes_arr=normalize_array(es_arr)
        n_p=n_p_value(es,es_arr)
        FDR=n_p_value(nes,nes_arr)

        print go_id," ",n_p
        gene_rankings.append([cluster, go_id, n_p,  FDR, es, nes])
    # prints out the rankings and significant values

    return gene_rankings

def normalize_score(es, es_arr):
    total = 0
    count=0
    for e in es_arr:
        if (e>0 and es>0) or (e<0 and es<0):
            total += e
            count+=1
    return es / (total/count)
def normalize_array(es_arr):

    null_distr=[]
    total = 0
    count = 0
    for es in es_arr:
        for e in es_arr:
            if (e>0 and es>0) or (e<0 and es<0):
                total += e
                count+=1
        null_distr.append(es / (total/count))
    return null_distr

def n_p_value(es, es_arr):
    total = 0
    for e in es_arr:
        total = total + e
    mean = total/len(es_arr)
    if es == mean:
        return 1.0
    elif es >= mean:
        tail = []
        for e in es_arr:
            if e >= es:
                tail.append(e)
        return float(len(tail))/float(len(es_arr))
    else:
        tail = []
        for e in es_arr:
            if e <= es:
                tail.append(e)
        return float(len(tail))/float(len(es_arr))



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

def enrichment_test(test_name, print_option, mat=None, anno=None,  rate=None, output=None, cluster=None):
    use_parsers = False
    if mat == None:
        p = Parsers("-a -c -o -l -r")
        anno = GMT(p.args.annotation_list)
        mat = MAT(p.args.cluster_list)
        cluster=p.args.cluster_number
        use_parsers = True
    else:
        mat = MAT(mat)
        anno = GMT(anno)

    if test_name=="wilcoxon":
        gene_rankings=wilcoxon(mat, anno)
    elif test_name=="page":
        gene_rankings=page(mat, anno)
    elif test_name=="gsea":
        gene_rankings=gsea(mat, anno, cluster)

    # prints out the rankings and significant values
    if use_parsers:
        return E_OUT(gene_rankings, p.args.output, p.args.rate, mat, anno).printout(print_option)
    else:
        return E_OUT(gene_rankings, output, rate, mat, anno).printout(print_option)


if __name__ == '__main__':

    #choose fisher_exact, chi_squared, hypergeometric, or binomial
    enrichment_test("gsea", True)