# class for generating enrichment test outputs
class OUT:
    def __init__(self, all_rankings, significant_rankings, output):
        self._all_rankings = all_rankings
        self._output = open(output, "r+")
        self._significant_rankings = significant_rankings

    # printing function for gsea
    # significant parameter is a boolean specifying if only significant results are desired according to the alpha level
    def printout_GSEA(self, print_to_console, significant_only):

        # print all significant gene sets

        if print_to_console:
            print("\ncluster\texpr_list.ngenes\tannotation_id\tannotation.ngenes\tp_value\tFDR\tes\tnes")
        self._output.write("\ncluster\texpr_list.ngenes\tannotation_id\tannotation.ngenes\tp_value\tFDR\tes\tnes")

        if significant_only:
            rankings = self._significant_rankings
        else:
            rankings = self._all_rankings

        for i, E_Result in enumerate(rankings):
            output_arr = []
            output_arr.append(E_Result.expr_cluster)
            output_arr.append(E_Result.expr_list_ngenes)
            output_arr.append(E_Result.anno_id)
            output_arr.append(E_Result.anno_ngenes)
            output_arr.append(E_Result.p_value)
            output_arr.append(E_Result.FDR)
            output_arr.append(E_Result.es)
            output_arr.append(E_Result.nes)
            output_arr = map(str, output_arr)

            self._output.write('\t'.join(output_arr) + "\n")

            if print_to_console:
                print('\t'.join(output_arr))

    # printing function for non-gsea tests
    # significant parameter is a boolean specifying if only significant results are desired according to the alpha level
    def printout_E(self, print_to_console, significant_only):

        # print all significant gene sets
        if print_to_console:
            print("\ncluster\texpr_list.ngenes\tannotation_id\tannotation.ngenes\tp_value\tFDR")
        self._output.write("\ncluster\texpr_list.ngenes\tannotation_id\tannotation.ngenes\tp_value\tFDR")

        if significant_only:
            rankings = self._significant_rankings
        else:
            rankings = self._all_rankings

        for i, E_Result in enumerate(rankings):
            output_arr = []
            output_arr.append(E_Result.expr_cluster)
            output_arr.append(E_Result.expr_list_ngenes)
            output_arr.append(E_Result.anno_id)
            output_arr.append(E_Result.anno_ngenes)
            output_arr.append(E_Result.p_value)
            output_arr.append(E_Result.FDR)
            output_arr = map(str, output_arr)

            self._output.write('\t'.join(output_arr) + "\n")

            if print_to_console:
                print('\t'.join(output_arr))
