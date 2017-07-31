'''
File name: overrep_unittests.py
Authors: Jimmy Shi
Date created: 6/28/2017
Date last modified:
Python Version: 2.7
'''

import os
import unittest
from multiprocessing import cpu_count

from elib.core.overrep_tests import *
from elib.utils.background import BACKGROUND


class TestStattests(unittest.TestCase):
    def test_gen_table(self):
        anno = GMT(os.path.join("files","unittest_files", "test_go.gmt")).genesets['0']
        sample = GMT(os.path.join("files","unittest_files", "test_gmt.gmt")).genesets['0']
        background = BACKGROUND([], os.path.join("files","unittest_files", "test_background.txt"))

        # contingency table from http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf
        # slide 24-26
        cont_table = gen_table(sample, anno, background)
        self.assertEqual(cont_table[0][0], 3)
        self.assertEqual(cont_table[0][1], 40)
        self.assertEqual(cont_table[1][0], 297)
        self.assertEqual(cont_table[1][1], 19960)

    def test_generate_inputs(self):
        anno = GMT(os.path.join("files","unittest_files", "test_go.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "test_background.txt"))
        self.assertEqual(generate_inputs(anno, background)[0][0], '0')
        test_set = list(range(0, 40))
        test_set = set([str(i) for i in test_set])
        input_items = generate_inputs(anno, background)
        self.assertEqual(input_items[0][1], test_set)
        self.assertEqual(input_items[0][2], background)

    def test_benjamini_hochberg(self):
        rankings = [OverrepResult(0, 0, 0, 0, 0.000162523456526, 0, 0), OverrepResult(0, 0, 0, 0, 0.0154958566105, 0, 0)
            , OverrepResult(0, 0, 0, 0, 0.021212455155, 0, 0), OverrepResult(0, 0, 0, 0, 0.0385370130677, 0, 0)
            , OverrepResult(0, 0, 0, 0, 0.118513682097, 0, 0), OverrepResult(0, 0, 0, 0, 0.138018032079, 0, 0)
            , OverrepResult(0, 0, 0, 0, 0.297867156785, 0, 0)]

        BH_result = benjamini_hochberg(rankings)
        self.assertAlmostEqual(BH_result[0].FDR, 0.001137664, delta=0.0001)
        self.assertAlmostEqual(BH_result[2].FDR, 0.04949573, delta=0.0001)
        self.assertAlmostEqual(BH_result[4].FDR, 0.1610210, delta=0.0001)

    def test_significance_filter(self):
        rankings = [OverrepResult(0, 0, 0, 0, 0, 0, 0.006), OverrepResult(1, 0, 0, 0, 0, 0, 0.03),
                    OverrepResult(2, 0, 0, 0, 0, 0, 0.2), OverrepResult(3, 0, 0, 0, 0, 0, 0.05)]
        rankings = significance_filter(rankings, 0.05)
        test_arr = []
        for i, item in enumerate(rankings):
            test_arr.append(item.gsid)
        self.assertTrue(0 in test_arr)
        self.assertTrue(1 in test_arr)
        self.assertTrue(2 not in test_arr)
        self.assertTrue(3 in test_arr)

    def test_binomial(self):
        anno = GMT(os.path.join("files","unittest_files", "GO_shortened.gmt"))
        sample = GMT(os.path.join("files","unittest_files", "GMT.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "BACKGROUND.txt"))

        test_result = binomial(sample, anno, 1, background, cpu_count())
        self.assertAlmostEqual(float(test_result[0][0].p_value), 0.000147087950101, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][0].FDR), 0.00367719875252, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][1].p_value), 0.0148019510684, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][1].FDR), 0.173734101344, delta=0.0001)

        anno = GMT(os.path.join("files","unittest_files", "test_go.gmt"))
        sample = GMT(os.path.join("files","unittest_files", "test_gmt.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "test_background.txt"))

        # results from http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf
        # slide 25
        self.assertAlmostEqual(float(binomial(sample, anno, 1, background, cpu_count())[0][0].p_value), 0.0229768421702,
                               delta=0.0001)

    def test_fisher(self):
        anno = GMT(os.path.join("files","unittest_files", "GO_shortened.gmt"))
        sample = GMT(os.path.join("files","unittest_files", "GMT.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "BACKGROUND.txt"))

        test_result = fisher_exact(sample, anno, 1, background, cpu_count())
        self.assertAlmostEqual(float(test_result[0][0].p_value), 0.000162523456526, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][0].FDR), 0.00406308641316, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][1].p_value), 0.0154958566105, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][1].FDR), 0.176770459625, delta=0.0001)

        anno = GMT(os.path.join("files","unittest_files", "test_go.gmt"))
        sample = GMT(os.path.join("files","unittest_files", "test_gmt.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "test_background.txt"))

        # results from http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf
        # slide 24
        self.assertAlmostEqual(float(fisher_exact(sample, anno, 1, background, cpu_count())[0][0].p_value), 0.0255246814673,
                               delta=0.0001)

    def test_chi_squared(self):
        anno = GMT(os.path.join("files","unittest_files", "GO_shortened.gmt"))
        sample = GMT(os.path.join("files","unittest_files", "GMT.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "BACKGROUND.txt"))

        test_result = chi_squared(sample, anno, 1, background, cpu_count())
        self.assertAlmostEqual(float(test_result[0][0].p_value), 2.463009591e-12, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][0].FDR), 6.15752397749e-11, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][1].p_value), 1.07544595815e-10, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][1].FDR), 8.96204965125e-10, delta=0.0001)

        anno = GMT(os.path.join("files","unittest_files", "test_go.gmt"))
        sample = GMT(os.path.join("files","unittest_files", "test_gmt.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "test_background.txt"))

        self.assertAlmostEqual(float(chi_squared(sample, anno, 1, background, cpu_count())[0][0].p_value), 1.27701446634e-64,
                               delta=0.0001)

    def test_hypergeometric(self):
        anno = GMT(os.path.join("files","unittest_files", "GO_shortened.gmt"))
        sample = GMT(os.path.join("files","unittest_files", "GMT.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "BACKGROUND.txt"))

        test_result = hypergeometric(sample, anno, 1, background, cpu_count())
        self.assertAlmostEqual(float(test_result[0][0].p_value), 0.000139835378986, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][0].FDR), 0.00349588447465, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][1].p_value), 0.0148068629756, delta=0.0001)
        self.assertAlmostEqual(float(test_result[0][1].FDR), 0.173762853283, delta=0.0001)

        anno = GMT(os.path.join("files","unittest_files", "test_go.gmt"))
        sample = GMT(os.path.join("files","unittest_files", "test_gmt.gmt"))
        background = BACKGROUND([], os.path.join("files","unittest_files", "test_background.txt"))

        # results from http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf
        # slide 26
        self.assertAlmostEqual(float(hypergeometric(sample, anno, 1, background, cpu_count())[0][0].p_value), 0.0219349067622,
                               delta=0.0001)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStattests)
    unittest.TextTestRunner(verbosity=2).run(suite)
