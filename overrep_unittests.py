import unittest
import os
from flib.core.gmt import GMT

from overrep_tests import *
from utilities.background import BACKGROUND


class TestStattests(unittest.TestCase):


    def test_gen_table(self):
        anno = GMT(os.path.join("unittest_files","test_go.gmt")).genesets['0']
        sample = GMT(os.path.join("unittest_files","test_gmt.gmt")).genesets['0']
        background = BACKGROUND([],os.path.join("unittest_files","test_background.txt"))

        self.assertEqual(gen_table(sample, anno, background)[0][0], 3)
        self.assertEqual(gen_table(sample, anno, background)[0][1], 40)
        self.assertEqual(gen_table(sample, anno, background)[1][0], 297)
        self.assertEqual(gen_table(sample, anno, background)[1][1], 19960)

    def test_generate_inputs(self):
        anno = GMT(os.path.join("unittest_files","test_go.gmt"))
        background = BACKGROUND([],os.path.join("unittest_files","test_background.txt"))
        self.assertEqual(generate_inputs(anno, background)[0][0], '0')
        test_set=list(range(0,40))
        test_set=set([str(i) for i in test_set])
        self.assertEqual(generate_inputs(anno, background)[0][1], test_set)
        self.assertEqual(generate_inputs(anno, background)[0][2], background)

    def test_benjamini_hochberg(self):
        rankings=[EnrichmentResult(0,0,0,0,0.000162523456526,0,0),EnrichmentResult(0,0,0,0,0.0154958566105,0,0)
            ,EnrichmentResult(0,0,0,0,0.021212455155,0,0),EnrichmentResult(0,0,0,0,0.0385370130677,0,0)
            ,EnrichmentResult(0,0,0,0,0.118513682097,0,0),EnrichmentResult(0,0,0,0,0.138018032079,0,0)
            ,EnrichmentResult(0,0,0,0,0.297867156785,0,0)]
        self.assertAlmostEqual(benjamini_hochberg(rankings)[0].FDR, 0.001137664, delta=0.0001)
        self.assertAlmostEqual(benjamini_hochberg(rankings)[2].FDR, 0.04949573, delta=0.0001)
        self.assertAlmostEqual(benjamini_hochberg(rankings)[4].FDR, 0.1610210, delta=0.0001)

    def test_significance_filter(self):
        rankings = [EnrichmentResult(0, 0, 0, 0, 0,0, 0.006), EnrichmentResult(1, 0, 0, 0,0,0,0.03), EnrichmentResult(2, 0, 0, 0,0,0, 0.2)]
        rankings=significance_filter(rankings,0.05)
        test_arr=[]
        for i, item in enumerate(rankings):
            test_arr.append(item.gsid)
        self.assertTrue(0 in test_arr)
        self.assertTrue(1 in test_arr)
        self.assertTrue(2 not in test_arr)

    def test_binomial(self):

        anno = GMT(os.path.join("unittest_files","test_go.gmt"))
        sample = GMT(os.path.join("unittest_files","test_gmt.gmt"))
        background = BACKGROUND([],os.path.join("unittest_files","test_background.txt"))

        self.assertAlmostEqual(float(binomial(sample, anno, 1, background)[0][0].p_value), 0.0229768421702,delta=0.0001)

    def test_fisher(self):
        anno = GMT(os.path.join("unittest_files","test_go.gmt"))
        sample =  GMT(os.path.join("unittest_files","test_gmt.gmt"))
        background =  BACKGROUND([],os.path.join("unittest_files","test_background.txt"))

        self.assertAlmostEqual(float(fisher_exact(sample, anno, 1, background)[0][0].p_value),0.0255246814673,delta=0.0001)


    def test_chi_squared(self):
        anno = GMT(os.path.join("unittest_files","test_go.gmt"))
        sample = GMT(os.path.join("unittest_files","test_gmt.gmt"))
        background =  BACKGROUND([],os.path.join("unittest_files","test_background.txt"))

        self.assertAlmostEqual(float(chi_squared(sample, anno, 1, background)[0][0].p_value), 1.27701446634e-64,delta=0.0001)


    def test_hypergeometric(self):
        anno = GMT(os.path.join("unittest_files","test_go.gmt"))
        sample = GMT(os.path.join("unittest_files","test_gmt.gmt"))
        background =  BACKGROUND([],os.path.join("unittest_files","test_background.txt"))

        self.assertAlmostEqual(float(hypergeometric(sample, anno, 1, background)[0][0].p_value), 0.0219349067622,delta=0.0001)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStattests)
    unittest.TextTestRunner(verbosity=2).run(suite)