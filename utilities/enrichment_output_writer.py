'''
File name: enrichment_output_writer.py
Authors: Jimmy Shi
Date created: 6/28/2017
Date last modified:
Python Version: 2.7
'''

import os
import webbrowser
from HTML import table

# class for generating enrichment test outputs
class OUT:
    def __init__(self, all_rankings, significant_rankings, output):
        self._all_rankings = all_rankings
        self._output = open(output, "r+")
        self.deleteContent(self._output)
        self._significant_rankings = significant_rankings

    def printout_GSEA(self, print_to_console, significant_only, precision):
        '''
        printing function for gsea

        :param bool print_to_console: Specifies whether or not to print to console
        :param bool significant_onky: Specifies whether or not to only output significant items
        :return: Nothing, will only write to the specified output file and to the console if specified
        '''

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

            if precision != -1:
                p_value=round(E_Result.p_value,precision)
                FDR = round(E_Result.FDR, precision)
                es=round(E_Result.es,precision)
                nes=round(E_Result.nes,precision)
            else:
                p_value = E_Result.p_value
                FDR = E_Result.FDR
                es = E_Result.es
                nes = E_Result.nes
            output_arr.append(p_value)
            output_arr.append(FDR)
            output_arr.append(es)
            output_arr.append(nes)
            output_arr = map(str, output_arr)

            self._output.write('\t'.join(output_arr) + "\n")

            if print_to_console:
                print('\t'.join(output_arr))
        self._output.close()

    def html_table_GSEA(self, significant_only,precision):
        '''
        generates and opens a html file for gsea in the default browser

        :param bool significant_only: Specifies whether or not to only output significant items
        :param bool precision: Specifies the decimal precision of the output, precision = 3 --> 3 decimal places, -1 will remove decimal restriction
        :return: Nothing, will only write to the table.html file and open the file in a browser
        '''

        html_output=open(os.path.join("utilities","table.html"), "r+")
        self.deleteContent(html_output)
        if significant_only:
            rankings = self._significant_rankings
        else:
            rankings = self._all_rankings

        output_arr = []
        for E_Result in rankings:
            if precision!=-1:
                p_value = round(E_Result.p_value, precision)
                FDR = round(E_Result.FDR, precision)
                es=round(E_Result.es, precision)
                nes=round(E_Result.nes, precision)
            else:
                p_value = E_Result.p_value
                FDR = E_Result.FDR
                es = E_Result.es
                nes = E_Result.nes

            next_row=[E_Result.expr_cluster,E_Result.expr_list_ngenes,E_Result.anno_id,E_Result.anno_ngenes,p_value, FDR, es,nes]
            next_row=map(str, next_row)
            output_arr.append(next_row)

        html_output.write(table(output_arr, header_row=["Expr Cluster", "Expr List Size", "Anno ID", "Anno Size","P Value", "FDR" , "ES", "NES"]))
        html_output.close()

        path = os.path.abspath(os.path.join("utilities","table.html"))
        url = "file://"+path
        webbrowser.open(url)

    def printout(self, print_to_console, significant_only, precision):
        '''
        printing function for all other enrichment tests

        :param bool print_to_console: Specifies whether or not to print to console
        :param bool significant_only: Specifies whether or not to only output significant items
        :return: Nothing, will only write to the specified output file and to the console if specified
        '''

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

            if precision!=-1:
                p_value = round(E_Result.p_value, precision)
                FDR = round(E_Result.FDR, precision)
            else:
                p_value = E_Result.p_value
                FDR = E_Result.FDR
            output_arr.append(p_value)
            output_arr.append(FDR)
            output_arr = map(str, output_arr)

            self._output.write('\t'.join(output_arr) + "\n")

            if print_to_console:
                print('\t'.join(output_arr))
        self._output.close()


    def html_table(self, significant_only, precision):
        '''
        generates and opens a html file in the default browser

        :param bool significant_only: Specifies whether or not to only output significant items
        :param bool precision: Specifies the decimal precision of the output, precision = 3 --> 3 decimal places, -1 will remove decimal restriction
        :return: Nothing, will only write to the table.html file and open the file in a browser
        '''

        html_output=open(os.path.join("utilities","table.html"), "r+")
        self.deleteContent(html_output)
        if significant_only:
            rankings = self._significant_rankings
        else:
            rankings = self._all_rankings

        output_arr = []
        for E_Result in rankings:
            if precision!=-1:
                p_value = round(E_Result.p_value, precision)
                FDR = round(E_Result.FDR, precision)
            else:
                p_value = E_Result.p_value
                FDR = E_Result.FDR
            next_row=[E_Result.expr_cluster,E_Result.expr_list_ngenes,E_Result.anno_id,E_Result.anno_ngenes,p_value,FDR]
            next_row = map(str, next_row)
            output_arr.append(next_row)

        html_output.write(table(output_arr, header_row=["Expr Cluster", "Expr List Size", "Anno ID", "Anno Size","P Value", "FDR"]))
        html_output.close()

        path = os.path.abspath(os.path.join("utilities","table.html"))
        url = "file://"+path
        webbrowser.open(url)

    def deleteContent(self, file):
        file.seek(0)
        file.truncate()