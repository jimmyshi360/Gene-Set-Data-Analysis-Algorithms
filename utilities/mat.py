from collections import OrderedDict
from operator import itemgetter
from collections import defaultdict

class MAT:
    def __init__(self, filename=None):
        self._dictionary={}
        self._dict = defaultdict()
        if filename:
            matfile = open(filename)
            for line in matfile:
                tok = line.strip().split('\t')

                if tok[0] != "human_entrez":
                    self._dictionary[tok[0]]= tok[1:]
                    self._dict[tok[0]]= tok[1:]
        self._ordered_dict = OrderedDict(sorted(self._dictionary.items(), key=itemgetter(1)))

    @property
    def dict(self):
        return self._dict

    @property
    def ordered_dict(self):
        return self._ordered_dict

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