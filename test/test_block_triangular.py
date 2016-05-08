import unittest
import random

from   sage.all import SR, matrix
import fuchsia
from   fuchsia import (block_triangular_form, transform)

class Test(unittest.TestCase):
    def assertIsTriangular(t, M1, M2, x, T, B):
        t.assertEqual(M2, transform(M1, x, T))

    def test_block_triangular_form_1(t):
        # Test with no transformation needed
        M = matrix([
          [1, 0, 0],
          [2, 3, 0],
          [4, 5, 6]
        ])
        MM, T, B = block_triangular_form(M)
        t.assertEqual(MM, M)
        t.assertEqual(T, matrix.identity(3))
        t.assertEqual(sorted(B), [(0, 1), (1, 1), (2, 1)])

    def test_block_triangular_form_2(t):
        # Test with no transformation possible
        M = matrix([
          [1, 2, 3, 4],
          [5, 6, 7, 8],
          [9, 1, 2, 3],
          [4, 5, 6, 7]
        ])
        MM, T, B = block_triangular_form(M)
        t.assertEqual(MM, M)
        t.assertEqual(T, matrix.identity(4))
        t.assertEqual(B, [(0, 4)])

    def test_block_triangular_form_3(t):
        m = matrix([
          [1, 0, 1, 0],
          [0, 1, 0, 1],
          [1, 0, 1, 0],
          [0, 0, 0, 1]
        ])
        mt, tt, b = block_triangular_form(m)
        x = SR.var("dummy")
        t.assertEqual(mt, transform(m, x, tt))
        t.assertEqual(matrix([
          [1, 1, 0, 0],
          [1, 1, 0, 0],
          [0, 0, 1, 0],
          [0, 0, 1, 1]
        ]), mt)

    def test_block_triangular_form_4(t):
        M = matrix([
          [1, 2, 3, 0, 0, 0],
          [4, 5, 6, 0, 0, 0],
          [7, 8, 9, 0, 0, 0],
          [2, 0, 0, 1, 2, 0],
          [0, 2, 0, 3, 4, 0],
          [0, 0, 2, 0, 0, 1]
        ])
        x = SR.var("dummy")
        T = matrix.identity(6)[random.sample(xrange(6), 6),:]
        M = transform(M, x, T)
        MM, T, B = block_triangular_form(M)
        t.assertEqual(MM, transform(M, x, T))
        t.assertEqual(sorted(s for o, s in B), [1, 2, 3])
        for o, s in B:
            for i in xrange(s):
                for j in xrange(s):
                    MM[o + i, o + j] = 0
        for i in xrange(6):
            for j in xrange(i):
                MM[i, j] = 0
        t.assertEqual(MM, matrix(6))
