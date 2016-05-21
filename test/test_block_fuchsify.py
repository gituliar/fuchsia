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

    def test_block_fuchsify_5(t):
        x, eps = SR.var("x eps")
        M = matrix([
            [eps/(x - 1), 0, 0],
            [(x**3 + x**2 - 3*x - 3)/(2*x**3 - 2*x**2 - x + 1),
                2*eps/(x - 3), 0],
            [-(x**3 + 3*x**2 - 3*x + 2)/(2*x**3 + x**2 - 2*x) + 2*(x**3 - x**2 - x + 1)/(2*x + 1),
                -3*(x**3 + x - 1)/(3*x**3 + x - 2),
                4*eps/(x - 5)]
        ])
        MM, T = block_fuchsify(M, x, eps)
        t.assertTrue((MM - transform(M, x, T).simplify_rational()).is_zero())
        pranks = singularities(MM, x).values()
        t.assertEqual(pranks, [0]*len(pranks))

    def test_block_fuchsify_6(t):
        x, eps = SR.var("x eps")
        M = matrix([
            [0,0],
            [1/(x**2-x+1)**2, 0]
        ])
        MM, T = block_fuchsify(M, x, eps)
        t.assertTrue((MM - transform(M, x, T).simplify_rational()).is_zero())
        pranks = singularities(MM, x).values()
        t.assertEqual(pranks, [0]*len(pranks))
