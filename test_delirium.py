#!/usr/bin/env python

import unittest
from random import randint

from delirium import *

def randpoly(x, maxrank=3):
    return sum(randint(-3, 3)*x**i for i in range(maxrank + 1))

def randrat(x, maxrank=3):
    return randpoly(x, maxrank)/randpoly(x, maxrank)

def randpolym(x, size, maxrank=3):
    return Matrix([
        [randpoly(x, maxrank) for j in range(size)]
        for i in range(size)
    ])

def randratm(x, size, maxrank=3):
    return Matrix([
        [randrat(x, maxrank) for j in range(size)]
        for i in range(size)
    ])

class Test(unittest.TestCase):
    def test_transform_1(t):
        # transform(M, x, I) == M
        x = sage.symbolic.ring.var("x")
        M = randpolym(x, 3)
        MM = transform(M, x, identity_matrix(M.nrows()))
        t.assertEqual(MM, M)

    def test_transform_2(t):
        # transform(transform(M, x, I), x, I^-1) == M
        x = sage.symbolic.ring.var("x")
        M = randpolym(x, 2)
        T = randpolym(x, 2)
        invT = T.inverse()
        M1 = transform(M, x, T)
        M2 = transform(M1, x, invT)
        t.assertEqual(M2.simplify_rational(), M)

    def test_balance_1(t):
        # balance(P, x1, x2, x)*balance(P, x2, x1, x) == I
        x = sage.symbolic.ring.var("x")
        P = Matrix([[1, 1], [0, 0]])
        x1 = randint(-10, 10)
        x2 = randint(20, 30)
        b1 = balance(P, x1, x2, x)
        b2 = balance(P, x2, x1, x)
        t.assertEqual((b1*b2).simplify_rational(), identity_matrix(P.nrows()))
        t.assertEqual((b2*b1).simplify_rational(), identity_matrix(P.nrows()))

    def test_balance_2(t):
        # balance(P, x1, oo, x)*balance(P, oo, x1, x) == I
        x = sage.symbolic.ring.var("x")
        P = Matrix([[1, 1], [0, 0]])
        x1 = randint(-10, 10)
        b1 = balance(P, x1, oo, x)
        b2 = balance(P, oo, x1, x)
        t.assertEqual((b1*b2).simplify_rational(), identity_matrix(P.nrows()))
        t.assertEqual((b2*b1).simplify_rational(), identity_matrix(P.nrows()))
    
if __name__ == "__main__":
    unittest.main(verbosity=2)
