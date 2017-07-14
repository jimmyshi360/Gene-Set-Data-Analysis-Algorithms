import unittest

from flib.core.gmt import GMT

from overrep_tests import *
from utilities.background import BACKGROUND


class TestStattests(unittest.TestCase):


    def test_gen_table(self):
        anno = GMT("unittest_files\go.gmt").genesets['0']
        sample = GMT("unittest_files\gmt.gmt").genesets['0']
        background = BACKGROUND("unittest_files\\background.txt")

        self.assertEqual(gen_table(sample, anno, background)[0][0], 3)
        self.assertEqual(gen_table(sample, anno, background)[0][1], 40)
        self.assertEqual(gen_table(sample, anno, background)[1][0], 297)
        self.assertEqual(gen_table(sample, anno, background)[1][1], 19960)

    def test_build_inputs(self):
        anno = GMT("unittest_files\go.gmt")
        background = BACKGROUND("unittest_files\\background.txt")
        self.assertEqual(build_inputs(anno, background)[0][0], '0')
        test_set=list(range(0,40))
        test_set=set([str(i) for i in test_set])
        self.assertEqual(build_inputs(anno, background)[0][1], test_set)
        self.assertEqual(build_inputs(anno, background)[0][2], background)

    def test_binomial(self):
        anno = "unittest_files\go.gmt"
        sample = "unittest_files\gmt.gmt"
        background = "unittest_files\\background.txt"
        output= "unittest_files\output.txt"
        self.assertEqual(over_rep_test("binomial", False, sample, anno, background, 0.05, output)[0][5], "0.0229768421702")

    def test_fisher(self):

        anno = "unittest_files\go.gmt"
        sample =  "unittest_files\gmt.gmt"
        background = "unittest_files\\background.txt"
        output = "unittest_files\output.txt"
        self.assertEqual(over_rep_test("fisher_exact", False,sample, anno, background, 0.05, output)[0][5],"0.0255246814673")


    def test_chi_squared(self):

        anno = "unittest_files\go.gmt"
        sample = "unittest_files\gmt.gmt"
        background = "unittest_files\\background.txt"
        output = "unittest_files\output.txt"
        self.assertEqual(over_rep_test("chi_squared", False,sample, anno, background, 0.05, output)[0][5], "1.27701446634e-64")


    def test_hypergeometric(self):

        anno = "unittest_files\go.gmt"
        sample = "unittest_files\gmt.gmt"
        background = "unittest_files\\background.txt"
        output = "unittest_files\output.txt"

        self.assertEqual(
            over_rep_test("hypergeometric", False,sample, anno, background, 0.05, output)[0][5], "0.0219349067622")

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStattests)
    unittest.TextTestRunner(verbosity=2).run(suite)