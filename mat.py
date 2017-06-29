from collections import defaultdict

class MAT:
    def __init__(self, filename=None):
        self._matrix = defaultdict(set)
        if filename:
            matfile = open(filename)
            for line in matfile:
                tok = line.strip().split('\t')

                if tok[0] != "human_entrez":
                    (gid, genes) = tok[0], tok[1:]

                    self._matrix[gid] = genes
    @property
    def matrix(self):
        return self._matrix
