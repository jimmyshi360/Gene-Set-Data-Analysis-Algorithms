
class BACKGROUND:
    def __init__(self, filename=None):
        self._background_genes=[]
        if filename:
            matfile = open(filename)
            for line in matfile:
                tok = line.strip().split('\t')
                self._background_genes.append(tok[0])
    @property
    def background_genes(self):
        return self._background_genes
