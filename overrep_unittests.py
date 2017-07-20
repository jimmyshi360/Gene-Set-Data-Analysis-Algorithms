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

    def test_build_inputs(self):
        anno = GMT(os.path.join("unittest_files","test_go.gmt"))
        background = BACKGROUND([],os.path.join("unittest_files","test_background.txt"))
        self.assertEqual(generate_inputs(anno, background)[0][0], '0')
        test_set=list(range(0,40))
        test_set=set([str(i) for i in test_set])
        self.assertEqual(generate_inputs(anno, background)[0][1], test_set)
        self.assertEqual(generate_inputs(anno, background)[0][2], background)

    def test_binomial(self):


        anno = os.path.join("unittest_files","test_go.gmt")
        sample = os.path.join("unittest_files","test_gmt.gmt")
        background = os.path.join("unittest_files","test_background.txt")
        output= os.path.join("unittest_files","test_output.txt")
        self.assertAlmostEqual(float(over_rep_test("binomial", False, sample, anno, background, 0.05, output)[0][5]), 0.0229768421702,delta=0.00001)

    def test_fisher(self):
        anno = os.path.join("unittest_files","test_go.gmt")
        sample =  os.path.join("unittest_files","test_gmt.gmt")
        background = os.path.join("unittest_files","test_background.txt")
        output = os.path.join("unittest_files","test_output.txt")
        self.assertAlmostEqual(float(over_rep_test("fisher_exact", False,sample, anno, background, 0.05, output)[0][5]),0.0255246814673,delta=0.00001)


    def test_chi_squared(self):
        anno = os.path.join("unittest_files","test_go.gmt")
        sample = os.path.join("unittest_files","test_gmt.gmt")
        background = os.path.join("unittest_files","test_background.txt")
        output = os.path.join("unittest_files","test_output.txt")
        self.assertAlmostEqual(float(over_rep_test("chi_squared", False,sample, anno, background, 0.05, output)[0][5]), 1.27701446634e-64,delta=0.00001)


    def test_hypergeometric(self):
        anno = os.path.join("unittest_files","test_go.gmt")
        sample = os.path.join("unittest_files","test_gmt.gmt")
        background = os.path.join("unittest_files","test_background.txt")
        output = os.path.join("unittest_files","test_output.txt")

        self.assertAlmostEqual(
        float(over_rep_test("hypergeometric", False,sample, anno, background, 0.05, output)[0][5]), 0.0219349067622,delta=0.00001)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStattests)
    unittest.TextTestRunner(verbosity=2).run(suite)