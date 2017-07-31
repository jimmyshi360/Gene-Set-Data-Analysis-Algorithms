'''
File name: mat.py
Authors: Jimmy Shi, Partha Rao
Date created: 6/28/2017
Date last modified:
Python Version: 2.7
'''

from collections import OrderedDict
from operator import itemgetter
from collections import defaultdict

class MAT:
    def __init__(self, filename=None):
        self._dictionary={}
        self._dict = defaultdict()
        self._labels=[]
        if filename:
            matfile = open(filename)
            for line in matfile:
                tok = line.strip().split('\t')

                if tok[0] != "human_entrez":
                    self._dictionary[tok[0]]= tok[1:]
                    self._dict[tok[0]]= tok[1:]
                else:
                    for i in range (0,len(tok)):
                        self._labels.append(tok[i])

        self._ordered_dict = OrderedDict(sorted(self._dictionary.items(), key=itemgetter(1)))

    @property
    def dict(self):
        return self._dict

    @property
    def ordered_dict(self):
        return self._ordered_dict

    @property
    def labels(self):
        return self._labels

    def sort(self,column):
        self._ordered_dict = OrderedDict(sorted(self._dictionary.items(), key=itemgetter(column + 1)))

    def scores(self,column):
        score_arr=[]
        for item in self._ordered_dict:
            score_arr.append(self._ordered_dict[item][column])
        return score_arr
    def ids(self):
        id_arr=[]
        for item in self._ordered_dict:
            id_arr.append(item)
        return id_arr

    #EXPERIMENTAL CONVERTERS, NOT TESTED
    def mat_to_rnk(self, output):
        output = open(output, "r+")
        output.write("\n")

        for i in self._dict.keys():
            output.write(str(i) + "\t" + str(self._dict[i][0]))
            output.write("\n")

    def mat_to_gct(self, output):
        output = open(output, "r+")

        output.write("#1.2\n")
        output.write(str(len(self._dict)) + "\t" + str(len(self._dict['1'])) + "\n")
        output.write("NAME\tDescription\t")

        for i in range(0, len(self._labels)):
            if i != 0:
                output.write(self._labels[i] + "\t")
        output.write("\n")

        for i in self._dict.keys():
            if len(self._dict[i]) != 0:
                output.write(str(i) + "\t")
                for j in self._dict[i]:
                    output.write(str(j) + "\t")
                output.write("\n")