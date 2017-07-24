# class for generating enrichment test outputs
class OUT:
    def __init__(self, gene_rankings, significant_rankings, output):
        self._gene_rankings = gene_rankings
        self._output = open(output, "r+")
        self._significant_rankings=significant_rankings

    # printing function for gsea
    # significant parameter is a boolean specifying if only significant results are desired according to the alpha level
    def printout_GSEA(self, print_to_console, significant):

        # print all significant gene sets

        if print_to_console:
            print("\nSIGNIFICANT VALUES")
            print("\ncluster\trankings.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR\tes\tnes")
        self._output.write("\ncluster\trankings.ngenes\tgs2\tgs2.ngenes\tFDR\tpvalue\tes\tnes")

        if significant:
            rankings=self._significant_rankings
        else:
            rankings =self._gene_rankings

        for i,E_Result in enumerate(rankings):
            output_arr = []
            output_arr.append(E_Result.cluster)
            output_arr.append(E_Result.cluster_ngenes)
            output_arr.append(E_Result.go_id)
            output_arr.append(E_Result.gs2_ngenes)
            output_arr.append(E_Result.p_value)
            output_arr.append(E_Result.FDR)
            output_arr.append(E_Result.es)
            output_arr.append(E_Result.nes)
            output_arr = map(str, output_arr)

            self._output.write('\t'.join(output_arr) + "\n")

            if print_to_console:
                print('\t'.join(output_arr))

    #printing function for non-gsea tests
    #significant parameter is a boolean specifying if only significant results are desired according to the alpha level
    def printout_E(self, print_to_console, significant):

        # print all significant gene sets
        if print_to_console:
            print("\nSIGNIFICANT VALUES")
            print("\ncluster\tcluster.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR")
        self._output.write("\ncluster\tcluster.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR")

        if significant:
            rankings=self._significant_rankings
        else:
            rankings =self._gene_rankings

        for i,E_Result in enumerate(rankings):
            output_arr = []
            output_arr.append(E_Result.cluster)
            output_arr.append(E_Result.cluster_ngenes)
            output_arr.append(E_Result.go_id)
            output_arr.append(E_Result.gs2_ngenes)
            output_arr.append(E_Result.p_value)
            output_arr.append(E_Result.FDR)
            output_arr = map(str, output_arr)

            self._output.write('\t'.join(output_arr) + "\n")

            if print_to_console:
                print('\t'.join(output_arr))



