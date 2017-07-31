import os
import unittest

from enrichment_tests import *
from flib.core.gmt import GMT
from multiprocessing import cpu_count

class TestStattests(unittest.TestCase):

    def test_generate_inputs(self):
        anno = GMT(os.path.join("unittest_files", "test_go.gmt"))
        expr_list = MAT(os.path.join("test_files", "CLUSTERS.mat"))
        cluster = 0
        permutations = 5

        test_set = list(range(0, 40))
        test_set = set([str(i) for i in test_set])

        inputs=generate_inputs(anno, cluster, expr_list, permutations)
        self.assertEqual(inputs[0].anno_id, '0')
        self.assertEqual(inputs[0].anno_list, test_set)
        self.assertEqual(inputs[0].expr_cluster, 0)
        self.assertEqual(inputs[0].expr_list, expr_list)
        self.assertEqual(inputs[0].permutations, 5)

    def test_gsea(self):
        anno = GMT(os.path.join("unittest_files", "GO_shortened.gmt"))
        expr_list = MAT(os.path.join("test_files", "CLUSTERS.mat"))
        cluster = 0
        permutations = 5

        gsea_result = gsea(expr_list, cluster, anno, permutations, 1.0, 1, cpu_count())
        self.assertAlmostEqual(float(gsea_result[0][0].es), 0.629079429538, delta=0.0001)
        self.assertAlmostEqual(float(gsea_result[0][1].es), 0.458915389493, delta=0.0001)

    def test_enrichment_score(self):
        anno = GMT(os.path.join("unittest_files", "GO_shortened.gmt")).genesets['GO:0070507']
        expr_list = MAT(os.path.join("test_files", "CLUSTERS.mat"))
        cluster = 0

        self.assertAlmostEqual(float(enrichment_score(anno, cluster, expr_list, 1)), 0.629079429538, delta=0.0001)

    def test_normalize_score(self):
        arr = [0.40078733453402793, 0.32961331889953255, 0.3281015607792101, 0.3156060188459115, 0.2882993457986646]
        self.assertAlmostEqual(normalize_score(0.629079429538, arr), 1.89207339265, delta=0.0001)

    def test_wilcoxon(self):
        anno = GMT(os.path.join("unittest_files", "GO_shortened.gmt"))
        expr_list = MAT(os.path.join("test_files", "CLUSTERS.mat"))
        cluster = 0

        wilcoxon_result=wilcoxon(expr_list, cluster, anno, 1, cpu_count())
        self.assertAlmostEqual(float(wilcoxon_result[0][0].p_value), 1.1206994619e-10, delta=0.0001)
        self.assertAlmostEqual(float(wilcoxon_result[0][1].p_value), 0.002463584445950147, delta=0.0001)

    def test_page(self):
        anno = GMT(os.path.join("unittest_files", "GO_shortened.gmt"))
        expr_list = MAT(os.path.join("test_files", "CLUSTERS.mat"))
        cluster = 0

        page_result=page(expr_list, cluster, anno, 1, cpu_count())
        self.assertAlmostEqual(float(page_result[0][0].p_value), 6.35801240607e-17, delta=0.0001)
        self.assertAlmostEqual(float(page_result[0][1].p_value), 1.1116067475e-10, delta=0.0001)

    def test_benjamini_hochberg(self):
        rankings = [EnrichmentResult(0, 0, 0, 0, 0.000162523456526, 0, 0),
                    EnrichmentResult(0, 0, 0, 0, 0.0154958566105, 0, 0)
            , EnrichmentResult(0, 0, 0, 0, 0.021212455155, 0, 0), EnrichmentResult(0, 0, 0, 0, 0.0385370130677, 0, 0)
            , EnrichmentResult(0, 0, 0, 0, 0.118513682097, 0, 0), EnrichmentResult(0, 0, 0, 0, 0.138018032079, 0, 0)
            , EnrichmentResult(0, 0, 0, 0, 0.297867156785, 0, 0)]

        BH_result=benjamini_hochberg(rankings)
        self.assertAlmostEqual(BH_result[0].FDR, 0.001137664, delta=0.0001)
        self.assertAlmostEqual(BH_result[2].FDR, 0.04949573, delta=0.0001)
        self.assertAlmostEqual(BH_result[4].FDR, 0.1610210, delta=0.0001)

    def test_significance_filter(self):
        rankings = [EnrichmentResult(0, 0, 0, 0, 0, 0.006), EnrichmentResult(0, 0, 1, 0, 0, 0.03),
                    EnrichmentResult(0, 0, 2, 0, 0, 0.2)]
        rankings = significance_filter(rankings, 0.05)
        test_arr = []
        for i, item in enumerate(rankings):
            test_arr.append(item.anno_id)
        self.assertTrue(0 in test_arr)
        self.assertTrue(1 in test_arr)
        self.assertTrue(2 not in test_arr)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStattests)
    unittest.TextTestRunner(verbosity=2).run(suite)
