import unittest

from   fuchsia.matrix import (
           alphabet, coefficient_low, degree_high, degree_low, import_matrix, matrix, var)

eps, x = var("eps, x")

class Test(unittest.TestCase):
    def test_alphabet(t):
        with open('examples/git_410.mtx') as f:
            m = import_matrix(f)
        t.assertEqual(
            alphabet(m, x),
            set([x-1, x, x+1]))

    def test_coefficient_low(t):
        with open('test/data/henn_324.mtx') as f:
            m = import_matrix(f)
        t.assertEqual(
            coefficient_low(m, 'x'),
            matrix([0, -1, 0, 0], 2, 2))

    def test_degree(t):
        # test 01
        with open('test/data/henn_324.mtx') as f:
            m = import_matrix(f)
        t.assertEqual(degree_high(m, 'x'), 0)
        t.assertEqual(degree_low(m, 'x'), -2)
        # test 02
        with open('test/data/bar_ex1.mtx') as f:
            m = import_matrix(f)
        t.assertEqual(degree_high(m, 'x'), 0)
        t.assertEqual(degree_low(m, 'x'), -9)

    def test_import_matrix(t):
        with open('test/data/henn_324.mtx') as f:
            m = import_matrix(f)
        t.assertEqual(m, matrix([eps/x, -1/x**2, 0, eps/(1+x)], 2, 2))

if __name__ == "__main__":
    unittest.main(verbosity=2)
