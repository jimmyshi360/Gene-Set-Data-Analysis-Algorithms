
# class for generating overrep test outputs
class OUT:
    def __init__(self, gene_rankings, significant_rankings, output, sample, anno):
        self._gene_rankings=gene_rankings
        self._output=open(output,"r+")
        self._sample=sample
        self._anno=anno
        self._significant_rankings=significant_rankings

    # printing function
    # significant parameter is a boolean specifying if only significant results are desired according to the alpha level
    def printout(self, print_option, significant):

        # print all significant gene sets
        if print_option:
            print("\nSIGNIFICANT VALUES")
            print("\ngs1\tgs1.ngenes\tgs2\tgs2.ngenes\tncommon\tpvalue\tFDR")
        self._output.write("\ngs1\tgs1.ngenes\tgs2\tgs2.ngenes\tncommon\tpvalue\tFDR\n")

        if significant:
            rankings=self._significant_rankings
        else:
            rankings=self._gene_rankings

        for i, row in enumerate(rankings):
            gsid=str(row[0])
            gs1_ngenes=str(len(self._sample.genesets[row[0]]))
            go_id=str(row[1])
            gs2_ngenes=str(len(self._anno.genesets[row[1]]))
            p_value=str(row[2])
            ncommon=str(row[3])
            FDR=str(row[4])
            self._output.write(gsid + "\t" + gs1_ngenes+"\t"+go_id+"\t"+gs2_ngenes+"\t"+ncommon+"\t"+p_value+"\t"+FDR+"\n")

            if print_option:
                print(gsid + "\t" + gs1_ngenes+"\t"+go_id+"\t"+gs2_ngenes+"\t"+ncommon+"\t"+p_value+"\t"+FDR)

