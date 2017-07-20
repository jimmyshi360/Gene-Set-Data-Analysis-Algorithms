import unittest

from enrichment_tests import *
import os
from flib.core.gmt import GMT

from utilities.background import BACKGROUND

class TestStattests(unittest.TestCase):


    def test_gsea(self):
        anno = GMT(os.path.join("unittest_files","GO_shortened.gmt"))
        rankings = MAT(os.path.join("test_files","CLUSTERS.mat"))
        cluster=0
        permutations=5
        self.assertAlmostEqual(float(gsea(rankings,anno,cluster,permutations)[0][4]), 0.629079429538, delta=0.0001)
        self.assertAlmostEqual(float(gsea(rankings,anno,cluster,permutations)[1][4]), 0.458915389493, delta=0.0001)

    def test_generate_inputs(self):
        anno = GMT(os.path.join("unittest_files", "test_go.gmt"))
        rankings = MAT(os.path.join("test_files", "CLUSTERS.mat"))
        cluster = 0
        permutations = 5

        test_set = list(range(0, 40))
        test_set = set([str(i) for i in test_set])

        self.assertEqual(generate_inputs(anno, cluster, rankings, permutations)[0][0],'0')
        self.assertEqual(generate_inputs(anno, cluster, rankings, permutations)[0][1], test_set)
        self.assertEqual(generate_inputs(anno, cluster, rankings, permutations)[0][2], 0)
        self.assertEqual(generate_inputs(anno, cluster, rankings, permutations)[0][3], rankings)
        self.assertEqual(generate_inputs(anno, cluster, rankings, permutations)[0][4], 5)
        self.assertEqual(len(generate_inputs(anno, cluster, rankings)[0]), 4)

    def test_enrichment_score(self):
        anno = GMT(os.path.join("unittest_files", "GO_shortened.gmt")).genesets['GO:0070507']
        rankings = MAT(os.path.join("test_files", "CLUSTERS.mat"))
        cluster = 0

        self.assertAlmostEqual(float(enrichment_score(anno,cluster,rankings,1)),  0.629079429538,delta=0.0001)

    def test_normalize_score(self):

        arr=[0.40078733453402793, 0.32961331889953255, 0.3281015607792101, 0.3156060188459115, 0.2882993457986646]
        self.assertAlmostEqual(normalize_score(0.629079429538,arr), 1.89207339265, delta=0.0001)

    def test_wilcoxon(self):
        anno = GMT(os.path.join("unittest_files","GO_shortened.gmt"))
        rankings = MAT(os.path.join("test_files","CLUSTERS.mat"))
        cluster=0
        self.assertAlmostEqual(float(wilcoxon(rankings,anno,cluster)[0][2]), 1.1206994619e-10, delta=0.0001)
        self.assertAlmostEqual(float(wilcoxon(rankings,anno,cluster)[1][2]), 0.0026701906744509285, delta=0.0001)

    def test_page(self):
        anno = GMT(os.path.join("unittest_files", "GO_shortened.gmt"))
        rankings = MAT(os.path.join("test_files", "CLUSTERS.mat"))
        cluster = 0
        self.assertAlmostEqual(float(page(rankings, anno, cluster)[0][2]), 6.35801240607e-17, delta=0.0001)
        self.assertAlmostEqual(float(page(rankings, anno, cluster)[1][2]), 1.1116067475e-10, delta=0.0001)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStattests)
    unittest.TextTestRunner(verbosity=2).run(suite)