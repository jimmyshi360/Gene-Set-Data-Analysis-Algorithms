
class OUT:
    def __init__(self, gene_rankings, output, alpha, sample, anno):
        self._gene_rankings=gene_rankings
        self._output=open(output,"r+")
        self._alpha=alpha
        self._sample=sample
        self._anno=anno
        self._significant_values=[]
        self._adjusted_p_values=[]

    def benjamini_hochberg(self):
        prev_bh_value = 0
        for i in range(0, len(self._gene_rankings)):
            bh_value = self._gene_rankings[i][2] * len(self._gene_rankings) / float(i + 1)
            bh_value = min(bh_value, 1)

            # to preserve monotonicity
            if prev_bh_value != 0:
                prev_bh_value = min(bh_value, prev_bh_value)
                if prev_bh_value< self._alpha:
                    self._significant_values.append(self._gene_rankings[i])
                    self._adjusted_p_values.append(prev_bh_value)
            if i == len(self._gene_rankings) - 1 and bh_value:
                self._significant_values.append(self._gene_rankings[i])
                self._adjusted_p_values.append(bh_value)
            prev_bh_value = bh_value

    def printout(self, print_option):
        #sorts the rankings in ascending order
        self._gene_rankings = sorted(self._gene_rankings, key=lambda line: float(line[0]))

        # grouping significant gene sets and multiple hypothesis test correction (hochberg)
        self.benjamini_hochberg()

        # print all significant gene sets
        if print_option:
            print("\nSIGNIFICANT VALUES")
            print("\ngs1\tgs1.ngenes\tgs2\tgs2.ngenes\tncommon\tpvalue\tFDR")
        self._output.write("\ngs1\tgs1.ngenes\tgs2\tgs2.ngenes\tncommon\tpvalue\tFDR\n")

        unit_test_arr=[]
        for i in range(0,len(self._significant_values)):
            x=self._significant_values[i]
            gsid=str(x[0])
            gs1_ngenes=str(len(self._sample.genesets[x[0]]))
            go_id=str(x[1])
            gs2_ngenes=str(len(self._anno.genesets[x[1]]))
            p_value=str(x[2])
            unit_test_arr.append([gsid,gs1_ngenes,go_id,gs2_ngenes,str(x[3]),p_value,str(self._adjusted_p_values[i])])
            self._output.write(gsid + "\t" + gs1_ngenes+"\t"+go_id+"\t"+gs2_ngenes+"\t"+str(x[3])+"\t"+p_value+"\t"+str(self._adjusted_p_values[i])+"\n")

            if print_option:
                print(gsid + "\t" + gs1_ngenes+"\t"+go_id+"\t"+gs2_ngenes+"\t"+str(x[3])+"\t"+p_value+"\t"+str(self._adjusted_p_values[i]))
        if print_option:
            print("LENGTH", len(unit_test_arr))
        return unit_test_arr
    def gen_out(self):
        self.printout()
        return self._gene_rankings