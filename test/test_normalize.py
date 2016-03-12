import unittest

from   sage.all import SR
from   fuchsia import is_normalized, matrix, matrix_c0, normalize, singularities, transform

class Test(unittest.TestCase):

    def test_is_normalized_1(t):
        x = SR.var("x")
        e = SR.var("epsilon")
        t.assertFalse(is_normalized(matrix([[1/x/2]]), x, e))
        t.assertFalse(is_normalized(matrix([[-1/x/2]]), x, e))
        t.assertTrue (is_normalized(matrix([[1/x/3]]), x, e))
        t.assertFalse(is_normalized(matrix([[x]]), x, e))
        t.assertFalse(is_normalized(matrix([[1/x**2]]), x, e))
        t.assertTrue (is_normalized(matrix([[(e+1/3)/x-1/2/(x-1)]]), x, e))

    def test_normalize_1(t):
        # Test with apparent singularities at 0 and oo, but not at 1.
        x = SR.var("x")
        M = matrix([
            [1/x, 5/(x-1), 0, 6/(x-1)],
            [0, 2/x, 0, 0],
            [0, 0, 3/x, 7/(x-1)],
            [6/(x-1), 0, 0, 1/x]
        ])

        N, T = normalize(M, x, SR.var("epsilon"))
        N = N.simplify_rational()
        t.assertEqual(N, transform(M, x, T).simplify_rational())
        for point, prank in singularities(N, x).iteritems():
            R = matrix_c0(N, x, point, prank)
            evlist = R.eigenvalues()
            t.assertEqual(evlist, [0]*len(evlist))

    def test_normalize_2(t):
        # Test with apparent singularities at 0, 1, and oo.
        x = SR.var("x")
        M = matrix([
            [1/x, 5/(x-1), 0, 6/(x-1)],
            [0, 2/(x-1), 0, 0],
            [0, 0, 3/x, 7/(x-1)],
            [6/(x-1), 0, 0, 1/x]
        ])

        N, T = normalize(M, x, SR.var("epsilon"))
        N = N.simplify_rational()
        t.assertEqual(N, transform(M, x, T).simplify_rational())
        for point, prank in singularities(N, x).iteritems():
            R = matrix_c0(N, x, point, prank)
            evlist = R.eigenvalues()
            t.assertEqual(evlist, [0]*len(evlist))

    def test_normalize_3(t):
        # Test with non-zero normalized eigenvalues
        x = SR.var("x")
        e = SR.var("epsilon")
        M = matrix([
            [(1-e)/x, 0],
            [0, (1+e)/3/x]
        ])

        N, T = normalize(M, x, e)
        N = N.simplify_rational()
        t.assertEqual(N, transform(M, x, T).simplify_rational())
        t.assertTrue(is_normalized(N, x, e))

    def test_normalize_4(t):
        # Test with non-zero normalized eigenvalues
        x, e = SR.var("x eps")
        M = matrix([
            [1/x/2, 0],
            [0, 0]
        ])

        with t.assertRaises(ValueError):
            N, T = normalize(M, x, e)

if __name__ == "__main__":
    unittest.main(verbosity=2)
