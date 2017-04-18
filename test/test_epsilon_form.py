import os.path
import unittest

from   sage.all import SR
from   fuchsia import (epsilon_form, import_matrix_from_file, is_fuchsian, is_normalized,
           transform, singularities)

class Test(unittest.TestCase):
    def assertReductionWorks(test, filename, fuchsian=False):
        m = import_matrix_from_file(filename)
        x, eps = SR.var("x eps")
        test.assertIn(x, m.variables())

        if not fuchsian:
            m_pranks = singularities(m, x).values()
            test.assertNotEqual(m_pranks, [0]*len(m_pranks))

        mt, t = epsilon_form(m, x, eps)
        test.assertTrue((mt-transform(m, x, t)).simplify_rational().is_zero())
        test.assertTrue(is_fuchsian(mt, x))
        test.assertTrue(is_normalized(mt, x, eps))
        test.assertNotIn(eps, (mt/eps).simplify_rational().variables())

    def test_git_409(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "data", "git_409.m"))

    def test_git_410(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "data", "git_410.m"))

    def test_henn_324(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "data", "henn_324.m"))

    def test_pap_3_50_slow(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "data", "pap_3_50.m"), fuchsian=True)
