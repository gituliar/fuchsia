import os.path
import unittest

from   sage.all import SR, matrix
from   fuchsia import (is_fuchsian, fuchsify_by_blocks, transform)

class Test(unittest.TestCase):
    def assert_fuchsify_by_blocks_works(test, m,b,x,eps):
        test.assertFalse(is_fuchsian(m, x))

        mt, t = fuchsify_by_blocks(m, b, x, eps)
        test.assertTrue((mt-transform(m, x, t)).simplify_rational().is_zero())
        test.assertTrue(is_fuchsian(mt, x))

    def test_fuchsify_by_blocks_01(test):
        x, eps = SR.var("x eps")
        m = matrix([
            [ eps/x,       0,       0],
            [     0, 2*eps/x,  -eps/x],
            [1/x**2,       0, 3*eps/x],
        ])
        b = [(0,1),(1,2)]
        test.assert_fuchsify_by_blocks_works(m, b, x, eps)

    def test_fuchsify_by_blocks_02(test):
        x, eps = SR.var("x eps")
        m = matrix([
            [ eps/x/(x-1),       0,       0],
            [  1/(x-1)**2, 2*eps/x,  -eps/x],
            [      1/x**2,       0, 3*eps/x],
        ])
        b = [(0,1),(1,2)]
        test.assert_fuchsify_by_blocks_works(m, b, x, eps)

    def test_fuchsify_by_blocks_03(test):
        x, eps = SR.var("x eps")
        m = matrix([
            [  eps/x,       0,       0],
            [ 1/x**2, 2*eps/x,       0],
            [ 1/x**2,  2/x**2, 3*eps/x],
        ])
        b = [(0,1),(1,1),(2,1)]
        test.assert_fuchsify_by_blocks_works(m, b, x, eps)

    def test_fuchsify_by_blocks_04(test):
        x, eps = SR.var("x eps")
        m = matrix([
            [     eps/x/(x-1),               0,             0],
            [ 1/x**2/(x-1)**2,   2*eps/x/(x-1),             0],
            [ 1/x**2/(x-1)**2, 2/x**2/(x-1)**2, 3*eps/x/(x-1)],
        ])
        b = [(0,1),(1,1),(2,1)]
        test.assert_fuchsify_by_blocks_works(m, b, x, eps)

    def test_fuchsify_by_blocks_05(test):
        x, eps = SR.var("x eps")
        m = matrix([
            [ eps/x,       0,       0],
            [     0, 2*eps/x,  -eps/x],
            [  x**2,       0, 3*eps/x],
        ])
        b = [(0,1),(1,2)]
        test.assert_fuchsify_by_blocks_works(m, b, x, eps)
