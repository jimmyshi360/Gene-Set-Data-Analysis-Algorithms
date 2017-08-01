'''
File name: overrep_output_writer.py
Authors: Jimmy Shi
Date created: 6/28/2017
Date last modified:
Python Version: 2.7
'''

import webbrowser
from HTML import table
import os

# class for generating overrep test outputs
class OUT:
    def __init__(self, all_rankings, significant_rankings, output, test_name):
        self._all_rankings = all_rankings
        self._output = open(output, "r+")
        self.deleteContent(self._output)
        self._significant_rankings = significant_rankings
        self._test_name=test_name

    def printout(self, print_to_console, significant_only, precision):
        '''
        printing function for overrrep tests

        :param bool print_to_console: Specifies whether or not to print to console
        :param bool significant_only: Specifies whether or not to only output significant items
        :return: Nothing, will only write to the specified output file and to the console if specified
        '''

        if print_to_console:
            print("\nrank\tgs1\tgs1.ngenes\tannotation_id\tannotation.ngenes\tncommon\tp_value\tFDR")
        self._output.write("\nrank\tgs1\tgs1.ngenes\tannotation_id\tannotation.ngenes\tncommon\tp_value\tFDR\n")

        if significant_only:
            rankings = self._significant_rankings
        else:
            rankings = self._all_rankings
        counter=1
        for OR_Result in rankings:

            output_arr = []
            output_arr.append(counter)
            output_arr.append(OR_Result.gsid)
            output_arr.append(OR_Result.sample_set_ngenes)
            output_arr.append(OR_Result.anno_id)
            output_arr.append(OR_Result.anno_ngenes)
            output_arr.append(OR_Result.overlaps)
            if precision!=-1:
                p_value=round(OR_Result.p_value, precision)
                FDR = round(OR_Result.FDR, precision)
            else:
                p_value = OR_Result.p_value
                FDR = OR_Result.FDR
            output_arr.append(p_value)
            output_arr.append(FDR)
            output_arr = map(str, output_arr)

            self._output.write('\t'.join(output_arr) + "\n")

            if print_to_console:
                print('\t'.join(output_arr))
            counter+=1
        self._output.close()

    def html_table(self, significant_only, precision):
        '''
        generates and opens a html file in the default browser

        :param bool significant_only: Specifies whether or not to only output significant items
        :param bool precision: Specifies the decimal precision of the output, precision = 3 --> 3 decimal places, -1 will remove decimal restriction
        :return: Nothing, will only write to the table.html file and open the file in a browser
        '''

        html_output=open(os.path.abspath(os.path.join("../","utils","table.html")), "r+")
        self.deleteContent(html_output)
        if significant_only:
            rankings = self._significant_rankings
        else:
            rankings = self._all_rankings

        output_arr = []
        counter=1
        for OR_Result in rankings:
            if precision!=-1:
                p_value=round(OR_Result.p_value, precision)
                FDR = round(OR_Result.FDR, precision)
            else:
                p_value = OR_Result.p_value
                FDR = OR_Result.FDR

            next_row=[counter,OR_Result.gsid,OR_Result.sample_set_ngenes,OR_Result.anno_id,OR_Result.anno_ngenes,OR_Result.overlaps,p_value,FDR]
            next_row= map(str,next_row)
            output_arr.append(next_row)
            counter+=1
        html_output.write("<p>" + self._test_name.upper() + "</p>")
        html_output.write(table(output_arr, header_row=[" ","GSID", "Set Size", "Anno ID", "Anno Size","Set Anno Overlaps","P Value", "FDR" ]))
        html_output.close()

        path = os.path.abspath(os.path.join("../","utils","table.html"))
        url = "file://"+path
        webbrowser.open(url)


    def deleteContent(self, file):
        file.seek(0)
        file.truncate()