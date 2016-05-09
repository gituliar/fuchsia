import os.path
import unittest

from   sage.all import SR
from   fuchsia import (block_triangular_form, import_matrix_from_file, is_normalized,
           is_normalized_by_blocks, normalize_by_blocks, transform, simplify_by_factorization,
           singularities)

class Test(unittest.TestCase):
    def assertNormalizeBlocksWorks(test, filename):
        x, eps = SR.var("x eps")

        m = import_matrix_from_file(filename)
        test.assertIn(x, m.variables())
        test.assertIn(eps, m.variables())
        test.assertFalse(is_normalized(m, x, eps))

        m_pranks = singularities(m, x).values()
        test.assertNotEqual(m_pranks, [0]*len(m_pranks))

        m, t, b = block_triangular_form(m)
        mt, tt = normalize_by_blocks(m, b, x, eps)
        t = t*tt
        test.assertTrue((mt-transform(m, x, t)).simplify_rational().is_zero())
        test.assertTrue(is_normalized_by_blocks(mt, b, x, eps))

    def test_git_409(test):
        test.assertNormalizeBlocksWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "git_409.mtx"))

    def test_git_410(test):
        test.assertNormalizeBlocksWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "git_410.mtx"))
