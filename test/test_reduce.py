import os.path
import unittest

from   sage.all import SR
from   fuchsia import (import_matrix_from_file, is_normalized, factorize,
           fuchsify, normalize, transform, simplify_by_factorization, singularities)

class Test(unittest.TestCase):
    def assertReductionWorks(t, filename):
        M = import_matrix_from_file(filename)
        x, eps = SR.var("x eps")
        t.assertIn(x, M.variables())
        M_pranks = singularities(M, x).values()
        t.assertNotEqual(M_pranks, [0]*len(M_pranks))

        #1 Fuchsify
        m, t1 = simplify_by_factorization(M, x)
        Mf, t2 = fuchsify(m, x)
        Tf = t1*t2
        t.assertTrue((Mf-transform(M, x, Tf)).simplify_rational().is_zero())
        Mf_pranks = singularities(Mf, x).values()
        t.assertEqual(Mf_pranks, [0]*len(Mf_pranks))

        #2 Normalize
        t.assertFalse(is_normalized(Mf, x, eps))
        m, t1 = simplify_by_factorization(Mf, x)
        for Mn, t2 in normalize(m, x, eps): pass
        Tn = t1*t2
        t.assertTrue((Mn-transform(Mf, x, Tn)).simplify_rational().is_zero())
        t.assertTrue(is_normalized(Mn, x, eps))

        #3 Factorize
        t.assertIn(eps, Mn.variables())
        m, t1 = simplify_by_factorization(Mn, x)
        Mc, t2 = factorize(m, x, eps, seed=3)
        Tc = t1*t2
        t.assertTrue((Mc-transform(Mn, x, Tc)).simplify_rational().is_zero())
        t.assertNotIn(eps, (Mc/eps).simplify_rational().variables())

    def test_git_409(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "git_409.mtx"))

    def test_git_410(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "git_410.mtx"))

    def test_henn_324(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "henn_324.mtx"))

if __name__ == "__main__":
    unittest.main(verbosity=2)
