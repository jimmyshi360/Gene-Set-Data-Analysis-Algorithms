


# class for generating overrep test outputs
class OUT:
    def __init__(self, gene_rankings, significant_rankings, output):
        self._gene_rankings=gene_rankings
        self._output=open(output,"r+")
        self._significant_rankings=significant_rankings

    # printing function
    # significant parameter is a boolean specifying if only significant results are desired according to the alpha level
    def printout(self, print_to_console, significant):

        # print all significant gene sets
        if print_to_console:
            print("\nSIGNIFICANT VALUES")
            print("\ngs1\tgs1.ngenes\tgs2\tgs2.ngenes\tncommon\tpvalue\tFDR")
        self._output.write("\ngs1\tgs1.ngenes\tgs2\tgs2.ngenes\tncommon\tpvalue\tFDR\n")

        if significant:
            rankings=self._significant_rankings
        else:
            rankings=self._gene_rankings

        for E_Result in rankings:
            output_arr=[]
            output_arr.append(E_Result.gsid)
            output_arr.append(E_Result.gs1_ngenes)
            output_arr.append(E_Result.go_id)
            output_arr.append(E_Result.gs2_ngenes)
            output_arr.append(E_Result.overlaps)
            output_arr.append(E_Result.p_value)
            output_arr.append(E_Result.FDR)
            output_arr=map(str,output_arr)

            self._output.write('\t'.join(output_arr)+"\n")

            if print_to_console:
                print('\t'.join(output_arr))

