# class for generating enrichment test outputs
class OUT:
    def __init__(self, gene_rankings, output, alpha, mat, anno):
        self._gene_rankings = gene_rankings
        self._output = open(output, "r+")
        self._alpha = alpha
        self._mat = mat
        self._anno = anno
        self._significant_values = []
        self._adjusted_p_values = []


    # FDR correction for multiple hypothesis testing
    def benjamini_hochberg(self):
        prev_bh_value = 0
        for i in range(0, len(self._gene_rankings)):
            bh_value = self._gene_rankings[i][2] * len(self._gene_rankings) / float(i + 1)
            bh_value = min(bh_value, 1)

            # to preserve monotonicity
            if prev_bh_value != 0:
                prev_bh_value = min(bh_value, prev_bh_value)
                if prev_bh_value < self._alpha:
                    self._significant_values.append(self._gene_rankings[i])
                    self._adjusted_p_values.append(prev_bh_value)
            if i == len(self._gene_rankings) - 1 and bh_value:
                self._significant_values.append(self._gene_rankings[i])
                self._adjusted_p_values.append(bh_value)
            prev_bh_value = bh_value

    # for GSEA FDR
    def gsea_FDR(self):

        for i in range(0, len(self._gene_rankings)):
            if self._gene_rankings[i][3] < self._alpha:
                self._significant_values.append(self._gene_rankings[i])

    # printing function
    def printout_GSEA(self, print_option):


        # grouping significant gene sets and multiple hypothesis test correction (hochberg)

        self.gsea_FDR()

        # print all significant gene sets
        if print_option:
            print("\nSIGNIFICANT VALUES")
            print("\nrankings\trankings.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR\tes\tnes")
        self._output.write("\nrankings\trankings.ngenes\tgs2\tgs2.ngenes\tFDR\tpvalue\tes\tnes")

        unit_test_arr = []
        for i in range(0, len(self._significant_values)):
            x = self._significant_values[i]
            rankings = str(x[0])
            rankings_ngenes = str(len(self._mat.matrix))
            go_id = str(x[1])
            gs2_ngenes = str(len(self._anno.genesets[x[1]]))
            p_value = str(x[2])
            unit_test_arr.append(
                [rankings, rankings_ngenes, go_id, gs2_ngenes, str(x[3]), p_value, str(self._significant_values[i])])
            self._output.write(rankings + "\t" + rankings_ngenes + "\t" + go_id + "\t" + gs2_ngenes + "\t" + str(
                x[3]) + "\t" + p_value + "\t" + str(self._significant_values[i]) + "\t" + str(x[4]) + "\t" + str(
                x[5]) + "\n")

            if print_option:
                print(rankings + "\t" + rankings_ngenes + "\t" + go_id + "\t" + gs2_ngenes + "\t" + str(
                    x[3]) + "\t" + p_value + "\t" + str(self._significant_values[i]) + "\t" + str(x[4]) + "\t" + str(
                    x[5]))
        if print_option:
            print("LENGTH", len(unit_test_arr))
        return unit_test_arr

    def printout_E(self, print_option):


        # grouping significant gene sets and multiple hypothesis test correction (hochberg)
        self.benjamini_hochberg()

        # print all significant gene sets
        if print_option:
            print("\nSIGNIFICANT VALUES")
            print("\ncluster\tcluster.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR")
        self._output.write("\ncluster\tcluster.ngenes\tgs2\tgs2.ngenes\tpvalue\tFDR")

        unit_test_arr = []
        for i in range(0, len(self._significant_values)):
            x = self._significant_values[i]
            cluster = str(x[0])
            cluster_ngenes = str(len(self._anno.genesets[x[0]]))
            go_id = str(x[1])
            gs2_ngenes = str(len(self._anno.genesets[x[1]]))
            p_value = str(x[2])
            FDR=str(self._adjusted_p_values[i])
            unit_test_arr.append(
                [cluster, cluster_ngenes, go_id, gs2_ngenes, p_value, FDR])
            self._output.write(cluster + "\t" + cluster_ngenes + "\t" + go_id + "\t" + gs2_ngenes + "\t" + p_value + "\t" + FDR + "\n")

            if print_option:
                print(cluster + "\t" + cluster_ngenes + "\t" + go_id + "\t" + gs2_ngenes + "\t" + p_value + "\t" + FDR)

        return unit_test_arr


    def gen_out(self):
        self.printout_GSEA()
        return self._gene_rankings
