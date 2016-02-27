#!/usr/bin/env python

import unittest

from fuchsia.reduction import *

class Test(unittest.TestCase):

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

if __name__ == "__main__":
    unittest.main(verbosity=2)
