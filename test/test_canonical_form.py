import os.path
import unittest

from   sage.all import SR
from   fuchsia import (canonical_form, import_matrix_from_file, is_fuchsian, is_normalized,
           transform, singularities)

class Test(unittest.TestCase):
    def assertReductionWorks(test, filename, fuchsian=False):
        m = import_matrix_from_file(filename)
        x, eps = SR.var("x eps")
        test.assertIn(x, m.variables())

        if not fuchsian:
            m_pranks = singularities(m, x).values()
            test.assertNotEqual(m_pranks, [0]*len(m_pranks))

        mt, t = canonical_form(m, x, eps)
        test.assertTrue((mt-transform(m, x, t)).simplify_rational().is_zero())
        test.assertTrue(is_fuchsian(mt, x))
        test.assertTrue(is_normalized(mt, x, eps))
        test.assertNotIn(eps, (mt/eps).simplify_rational().variables())

    def test_git_409(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "git_409.mtx"))

    def test_git_410(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "git_410.mtx"))

    def test_henn_324(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "henn_324.mtx"))

    def test_pap_03_52_slow(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "pap_03_52.mtx"), fuchsian=True)
