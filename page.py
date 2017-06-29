from scipy import stats
from flib.core.gmt import GMT
from mat import MAT
import numpy as np
import math

def page(gmt, mat, output,  cluster,false_discovery_rate):
    gene_rankings = []
    gene_mean = 0
    gene_sd = 0
    geneset_mean = 0
    geneset_size = 0

    score_arr = []

    for gsid in gmt._genesets:
        for gene in gmt._genesets[gsid]:
            row_arr=list(mat._matrix[gene])
            score_arr.append(row_arr[cluster+1])

    score_arr=np.array(score_arr).astype(np.float)
    gene_mean = np.mean(score_arr)
    gene_sd = np.std(score_arr)
    for gsid in gmt._genesets:
        geneset_size = len(gmt._genesets[gsid])
        z_score = (geneset_mean - gene_mean) * math.sqrt(geneset_size) / gene_sd
        p_value = stats.norm.sf(abs(z_score))
        gene_rankings.append([p_value, gsid])
    gene_rankings = sorted(gene_rankings, key=lambda line: float(line[0]))

    # grouping significant gene sets and multiple hypothesis test correction (hochberg)
    significant_values = []
    for i in range(0, len(gene_rankings)):
        if gene_rankings[i][0] < float(i) / len(gene_rankings) * false_discovery_rate:
            significant_values.append(gene_rankings[i])
    # print rankings and write to output file, reverses ascending array before sorting
    print("\n\nRANKINGS")
    for set_arr in gene_rankings[::-1]:
        output.write(set_arr[1] + ": " + str(set_arr[0]))
        print(set_arr[1] + ": " + str(set_arr[0]))

    # print all significant gene sets
    print("\n\nSIGNIFICANT VALUES")
    for x in significant_values:
        print(x[1] + " " + str(x[0]))

gmt = GMT("C:\\Users\\Jimmy\\Documents\\GENOMICS\\list.txt")
mat = MAT("C:\\Users\\Jimmy\\Documents\\GENOMICS\\clusters.mat")
output = open("C:\\Users\\Jimmy\\Documents\\GENOMICS\\output.txt", "r+")

page(gmt,mat, output, 0,0.05)