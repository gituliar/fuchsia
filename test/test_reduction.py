#!/usr/bin/env python

import unittest
from random import randint

from delirium.reduction import *

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
        x = SR.var("x")
        M = randpolym(x, 3)
        MM = transform(M, x, identity_matrix(M.nrows()))
        t.assertEqual(MM, M)

    def test_transform_2(t):
        # transform(transform(M, x, I), x, I^-1) == M
        x = SR.var("x")
        M = randpolym(x, 2)
        T = randpolym(x, 2)
        invT = T.inverse()
        M1 = transform(M, x, T)
        M2 = transform(M1, x, invT)
        t.assertEqual(M2.simplify_rational(), M)

    def test_balance_1(t):
        # balance(P, x1, x2, x)*balance(P, x2, x1, x) == I
        x = SR.var("x")
        P = Matrix([[1, 1], [0, 0]])
        x1 = randint(-10, 10)
        x2 = randint(20, 30)
        b1 = balance(P, x1, x2, x)
        b2 = balance(P, x2, x1, x)
        t.assertEqual((b1*b2).simplify_rational(), identity_matrix(P.nrows()))
        t.assertEqual((b2*b1).simplify_rational(), identity_matrix(P.nrows()))

    def test_balance_2(t):
        # balance(P, x1, oo, x)*balance(P, oo, x1, x) == I
        x = SR.var("x")
        P = Matrix([[1, 1], [0, 0]])
        x1 = randint(-10, 10)
        b1 = balance(P, x1, oo, x)
        b2 = balance(P, oo, x1, x)
        t.assertEqual((b1*b2).simplify_rational(), identity_matrix(P.nrows()))
        t.assertEqual((b2*b1).simplify_rational(), identity_matrix(P.nrows()))

    def test_reduce_at_one_point_1(t):
        x = SR.var("x")
        M0 = matrix([
            [1/x, 4, 0, 5],
            [0, 2/x, 0, 0],
            [0, 0, 3/x, 6],
            [0, 0, 0, 4/x]
        ])

        u = matrix([
            [0, Rational((3, 5)), Rational((4, 5)), 0],
            [Rational((5, 13)), 0, 0, Rational((12, 13))]
        ])
        M1 = transform(M0, x, balance(u.transpose()*u, 0, 1, x))
        M1 = M1.simplify_rational()

        u = matrix([[8, 0, 15, 0]])/17
        M2 = transform(M1, x, balance(u.transpose()*u, 0, 2, x))
        M2 = M2.simplify_rational()

        M2_sing = singularities(M2, x)
        t.assertIn(0, M2_sing)
        t.assertEqual(M2_sing[0], 3)

        M3, T23 = reduce_at_one_point(M2, x, 0, 3)
        M3 = M3.simplify_rational()
        t.assertEqual(M3, transform(M2, x, T23).simplify_rational())

        M3_sing = singularities(M3, x)
        t.assertIn(0, M3_sing)
        t.assertEqual(M3_sing[0], 2)

        M4, T34 = reduce_at_one_point(M3, x, 0, 2)
        M4 = M4.simplify_rational()
        t.assertEqual(M4, transform(M3, x, T34).simplify_rational())

        M4_sing = singularities(M4, x)
        t.assertIn(0, M4_sing)
        t.assertEqual(M4_sing[0], 1)
    
if __name__ == "__main__":
    unittest.main(verbosity=2)
