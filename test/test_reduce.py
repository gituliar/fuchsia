import os.path
import unittest

from   sage.all import SR
from   fuchsia import (import_matrix_from_file, is_normalized, factorize,
           fuchsify, normalize, partial_fraction, transform, singularities)

class Test(unittest.TestCase):
    def assertReductionWorks(t, filename):
        M = import_matrix_from_file(filename)
        x, eps = SR.var("x eps")
        t.assertIn(x, M.variables())
        M_pranks = singularities(M, x).values()
        t.assertNotEqual(M_pranks, [0]*len(M_pranks))

        #1 Fuchsify
        Mf, Tf = fuchsify(M, x)
        Mf = Mf.simplify_rational()
        t.assertEqual(Mf, transform(M, x, Tf).simplify_rational())
        Mf_pranks = singularities(Mf, x).values()
        t.assertEqual(Mf_pranks, [0]*len(Mf_pranks))

        #2 Normalize
        t.assertFalse(is_normalized(Mf, x, eps))
        Mn, Tn = normalize(Mf, x, eps)
        Mn = Mn.simplify_rational()
        t.assertEqual(Mn, transform(Mf, x, Tn).simplify_rational())
        t.assertTrue(is_normalized(Mn, x, eps))

        #3 Factorize
        Mn = partial_fraction(Mn,x)
        t.assertIn(eps, Mn.variables())
        Mc, Tc = factorize(Mn, x, eps, seed=3)
        Mc = Mc.simplify_rational()
        t.assertEqual(Mc, transform(Mn, x, Tc).simplify_rational())
        t.assertNotIn(eps, (Mc/eps).simplify_rational().variables())

    def test_git_409(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "git_409.mtx"))

    @unittest.skip("Part #3 of the test is hanging")
    def test_git_410(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "git_410.mtx"))

    def test_henn_324(t):
        t.assertReductionWorks(os.path.join(os.path.dirname(__file__),
            "..", "examples", "henn_324.mtx"))

if __name__ == "__main__":
    unittest.main(verbosity=2)
