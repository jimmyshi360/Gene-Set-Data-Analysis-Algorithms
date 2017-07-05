from scipy import stats
from flib.core.gmt import GMT
from output_generator import OUT
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

def enrichment_score(gene_set, cluster, mat, weight):

    N = len(gene_set)

    set_score_sum=0
    for id in gene_set:
        if len(list(mat[id]))!=0:
            set_score_sum+=(float(list(mat[id])[cluster]))**weight

    max_ES = -sys.maxint
    min_ES = sys.maxint
    net_sum = 0
    sum_arr = []
    for i in mat.keys():
        if i in gene_set:
            if len(list(mat[id])) != 0:
                net_sum += (float(list(mat[i])[cluster]))**weight/set_score_sum
        else:
            net_sum += -1.0/(len(mat)-len(gene_set))
        if net_sum > max_ES:
            max_ES = net_sum
        if net_sum < min_ES:
            min_ES = net_sum
        sum_arr.append(net_sum)

    '''
    line_arr=np.repeat(0,25000)
    plt.plot(sum_arr)
    plt.plot(line_arr)
    plt.ylabel("ES Plot")
    plt.show()
    '''
    return max_ES

def permute(mat, gene_set, cluster):
    mat_keys = np.random.permutation(list(mat.matrix.keys()))
    return mat_keys[0:len(gene_set)]

def es_distr(set,mat,cluster):
    es_arr=[]
    for i in range(0,100):
        es_arr.append(enrichment_score(permute(mat,set,cluster),cluster,mat.matrix,1))
    return es_arr

def gsea(sample, mat,cluster):

    gene_rankings=[]
    for gsid in sample.genesets:
        es=enrichment_score(sample.genesets[gsid], cluster, mat.matrix, 1)
        es_arr=es_distr(sample.genesets[gsid], mat, cluster)
        es=normalize_score(es,es_arr)
        es_arr=normalize_array(es_arr)
        n_p=n_p_value(es,es_arr)

        '''
        plt.plot(es_arr)
        plt.ylabel("ES Distribution Plot")
        plt.show()
        '''

        gene_rankings.append([n_p,gsid])
    # prints out the rankings and significant values
    OUT(gene_rankings,"OUTPUT.txt",0.05).printout()

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

gsea(GMT("GMT.gmt"), MAT("CLUSTERS.mat"), 0)