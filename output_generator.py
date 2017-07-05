
class OUT:
    def __init__(self, gene_rankings, output, false_discovery_rate):
        self._gene_rankings=gene_rankings
        self._output=open(output,"r+")
        self._false_discovery_rate=false_discovery_rate

    def printout(self):

        #sorts the rankings in ascending order
        self._gene_rankings = sorted(self._gene_rankings, key=lambda line: float(line[0]))

        # grouping significant gene sets and multiple hypothesis test correction (hochberg)
        significant_values = []
        for i in range(0, len(self._gene_rankings)):
            if  self._gene_rankings[i][0] < float(i) / len( self._gene_rankings) *  self._false_discovery_rate:
                significant_values.append( self._gene_rankings[i])

        # print rankings and write to output file, reverses ascending array before sorting
        self._output.write("\nRANKINGS")
        print("\nRANKINGS")
        for set_arr in  self._gene_rankings:
            self._output.write(str(set_arr[1]) + ": " + str(set_arr[0]))
            print(str(set_arr[1]) + ": " + str(set_arr[0]))

        significant_values = sorted(significant_values, key=lambda line: float(line[0]))
        # print all significant gene sets
        self._output.write("\nSIGNIFICANT VALUES")
        print("\nSIGNIFICANT VALUES")
        for x in significant_values:
            self._output.write(str(x[1]) + " " + str(x[0]))
            print(str(x[1]) + " " + str(x[0]))