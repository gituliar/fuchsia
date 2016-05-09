import unittest
from   sage.all import SR, matrix
from   fuchsia import \
        block_fuchsify, import_matrix_from_file, singularities, \
        transform

class Test(unittest.TestCase):
    def test_block_fuchsify_1(t):
        x, eps = SR.var("x eps")
        M = matrix([
            [eps/x, 0],
            [1/x/x, eps/x]
        ])
        MM, T = block_fuchsify(M, x, eps)
        t.assertEqual(MM, transform(M, x, T).simplify_rational())
        pranks = singularities(MM, x).values()
        t.assertEqual(pranks, [0]*len(pranks))

    def test_block_fuchsify_2(t):
        x, eps = SR.var("x eps")
        M = matrix([
            [eps/x, 0],
            [1, eps/x]
        ])
        MM, T = block_fuchsify(M, x, eps)
        t.assertEqual(MM, transform(M, x, T).simplify_rational())
        pranks = singularities(MM, x).values()
        t.assertEqual(pranks, [0]*len(pranks))

    def test_block_fuchsify_3(t):
        x, eps = SR.var("x eps")
        M = matrix([
            [eps/x, 0, 0],
            [1/x, eps/x, 0],
            [1/x**3, x**3, eps/x]
        ])
        MM, T = block_fuchsify(M, x, eps)
        t.assertEqual(MM, transform(M, x, T).simplify_rational())
        pranks = singularities(MM, x).values()
        t.assertEqual(pranks, [0]*len(pranks))

    def test_block_fuchsify_4(t):
        x, eps = SR.var("x eps")
        M = import_matrix_from_file("examples/git_410.mtx")
        MM, T = block_fuchsify(M, x, eps)
        t.assertEqual(MM, transform(M, x, T).simplify_rational())
        pranks = singularities(MM, x).values()
        t.assertEqual(pranks, [0]*len(pranks))
