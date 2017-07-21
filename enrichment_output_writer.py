# class for generating enrichment test outputs
class OUT:
    def __init__(self, gene_rankings, significant_rankings, output, mat, anno):
        self._gene_rankings = gene_rankings
        self._output = open(output, "r+")
        self._mat = mat
        self._anno = anno
        self._significant_rankings=significant_rankings

    # printing function for gsea
    # significant parameter is a boolean specifying if only significant results are desired according to the alpha level
    def printout_GSEA(self, print_option, significant):

        # print all significant gene sets

        if print_option:
            print("\nSIGNIFICANT VALUES")
            print("\ncluster\trankings.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR\tes\tnes")
        self._output.write("\ncluster\trankings.ngenes\tgs2\tgs2.ngenes\tFDR\tpvalue\tes\tnes")

        if significant:
            rankings=self._significant_rankings
        else:
            rankings =self._gene_rankings

        for i,row in enumerate(rankings):
            cluster = str(row[0])
            rankings_ngenes = str(len(self._mat.dict))
            go_id = str(row[1])
            gs2_ngenes = str(len(self._anno.genesets[row[1]]))
            p_value = str(row[2])
            FDR=str(row[3])
            es=str(row[4])
            nes=str(row[5])
            self._output.write(cluster + "\t" + rankings_ngenes + "\t" + go_id + "\t" + gs2_ngenes + "\t" + p_value + "\t" + FDR + "\t" + es + "\t" +nes+ "\n")

            if print_option:
                print(cluster + "\t" + rankings_ngenes + "\t" + go_id + "\t" + gs2_ngenes + "\t" + p_value + "\t" + FDR + "\t" + es + "\t" + nes)

    #printing function for non-gsea tests
    #significant parameter is a boolean specifying if only significant results are desired according to the alpha level
    def printout_E(self, print_option, significant):

        # print all significant gene sets
        if print_option:
            print("\nSIGNIFICANT VALUES")
            print("\ncluster\tcluster.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR")
        self._output.write("\ncluster\tcluster.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR")

        if significant:
            rankings=self._significant_rankings
        else:
            rankings =self._gene_rankings
        for i,row in enumerate(rankings):

            cluster = str(row[0])
            cluster_ngenes = str(len(self._anno.genesets[row[0]]))
            go_id = str(row[1])
            gs2_ngenes = str(len(self._anno.genesets[row[1]]))
            p_value = str(row[2])
            FDR=str(row[3])
            self._output.write(cluster + "\t" + cluster_ngenes + "\t" + go_id + "\t" + gs2_ngenes + "\t" + p_value + "\t" + FDR + "\n")

            if print_option:
                print(cluster + "\t" + cluster_ngenes + "\t" + go_id + "\t" + gs2_ngenes + "\t" + p_value + "\t" + FDR)



