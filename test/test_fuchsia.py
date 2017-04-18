import unittest
from   random import randint
from   StringIO import StringIO

from   sage.all import SR
import fuchsia
from   fuchsia import \
    balance, balance_transform, identity_matrix, import_matrix_mathematica, \
    import_matrix_matrixmarket, import_matrix_from_file, \
    export_matrix_mathematica, export_matrix_matrixmarket, factorize, \
    fuchsify, limit_fixed, matrix, matrix_complexity, oo, Rational, \
    reduce_at_one_point, simplify_by_jordanification, singularities, \
    transform, FuchsiaError

def randpoly(x, maxrank=3):
    return sum(randint(-3, 3)*x**i for i in range(maxrank + 1))

def randrat(x, maxrank=3):
    return randpoly(x, maxrank)/randpoly(x, maxrank)

def randpolym(x, size, maxrank=3):
    return matrix([
        [randpoly(x, maxrank) for j in range(size)]
        for i in range(size)
    ])

def randratm(x, size, maxrank=3):
    return matrix([
        [randrat(x, maxrank) for j in range(size)]
        for i in range(size)
    ])

class Test(unittest.TestCase):
    def test_import_matrix_from_file_1(t):
        x, eps = SR.var("x eps")
        m = import_matrix_from_file("test/data/henn_324.m")
        t.assertEqual(set(m.variables()), set([x, eps]))
        t.assertEqual(m, matrix([[eps/x, 0], [-1/x**2, eps/(x + 1)]]))

    def test_import_export_matrixmarket(t):
        a, b = SR.var("v1 v2")
        M = matrix([[1, a, b], [a + b, Rational((2, 3)), a/b]])
        fout = StringIO()
        export_matrix_matrixmarket(fout, M)
        MM = import_matrix_matrixmarket(StringIO(fout.getvalue()))
        t.assertEqual(M, MM)

    def test_import_export_mathematica(t):
        a, b = SR.var("v1 v2")
        M = matrix([[1, a, b], [a + b, Rational((2, 3)), a/b]])
        fout = StringIO()
        export_matrix_mathematica(fout, M)
        MM = import_matrix_mathematica(StringIO(fout.getvalue()))
        t.assertEqual(M, MM)

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
        P = matrix([[1, 1], [0, 0]])
        x1 = randint(-10, 10)
        x2 = randint(20, 30)
        b1 = balance(P, x1, x2, x)
        b2 = balance(P, x2, x1, x)
        t.assertEqual((b1*b2).simplify_rational(), identity_matrix(P.nrows()))
        t.assertEqual((b2*b1).simplify_rational(), identity_matrix(P.nrows()))

    def test_balance_2(t):
        # balance(P, x1, oo, x)*balance(P, oo, x1, x) == I
        x = SR.var("x")
        P = matrix([[1, 1], [0, 0]])
        x1 = randint(-10, 10)
        b1 = balance(P, x1, oo, x)
        b2 = balance(P, oo, x1, x)
        t.assertEqual((b1*b2).simplify_rational(), identity_matrix(P.nrows()))
        t.assertEqual((b2*b1).simplify_rational(), identity_matrix(P.nrows()))

    def test_balance_transform_1(t):
        x = SR.var("x")
        M = randpolym(x, 2)
        P = matrix([[1, 1], [0, 0]])
        x1 = randint(-10, 10)
        x2 = randint(20, 30)
        b1 = balance(P, x1, x2, x)

        M1 = balance_transform(M, P, x1, x2, x)
        M2 = transform(M, x, balance(P, x1, x2, x))
        t.assertEqual(M1.simplify_rational(), M2.simplify_rational())

        M1 = balance_transform(M, P, x1, oo, x)
        M2 = transform(M, x, balance(P, x1, oo, x))
        t.assertEqual(M1.simplify_rational(), M2.simplify_rational())

        M1 = balance_transform(M, P, oo, x2, x)
        M2 = transform(M, x, balance(P, oo, x2, x))
        t.assertEqual(M1.simplify_rational(), M2.simplify_rational())

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
        t.assertEqual(M2_sing[0], 2)

        M3, T23 = reduce_at_one_point(M2, x, 0, 2)
        M3 = M3.simplify_rational()
        t.assertEqual(M3, transform(M2, x, T23).simplify_rational())

        M3_sing = singularities(M3, x)
        t.assertIn(0, M3_sing)
        t.assertEqual(M3_sing[0], 1)

        M4, T34 = reduce_at_one_point(M3, x, 0, 1)
        M4 = M4.simplify_rational()
        t.assertEqual(M4, transform(M3, x, T34).simplify_rational())

        M4_sing = singularities(M4, x)
        t.assertIn(0, M4_sing)
        t.assertEqual(M4_sing[0], 0)

    def test_fuchsify_1(t):
        x = SR.var("x")
        M = matrix([
            [1/x, 5, 0, 6],
            [0, 2/x, 0, 0],
            [0, 0, 3/x, 7],
            [0, 0, 0, 4/x]
        ])

        u = matrix([
            [0, Rational((3, 5)), Rational((4, 5)), 0],
            [Rational((5, 13)), 0, 0, Rational((12, 13))]
        ])
        M = transform(M, x, balance(u.transpose()*u, 0, 1, x))
        M = M.simplify_rational()

        u = matrix([[8, 0, 15, 0]])/17
        M = transform(M, x, balance(u.transpose()*u, 0, 2, x))
        M = M.simplify_rational()

        Mx, T = fuchsify(M, x)
        Mx = Mx.simplify_rational()
        t.assertEqual(Mx, transform(M, x, T).simplify_rational())

        pranks = singularities(Mx, x).values()
        t.assertEqual(pranks, [0]*len(pranks))

    def test_fuchsify_2(t):
        x = SR.var("x")
        M = matrix([
            [0, 1/x/(x-1), 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ])
        u = matrix([[6, 3, 2, 0]])/7
        P = u.transpose()*u
        M = balance_transform(M, P, 1, 0, x).simplify_rational()
        M = balance_transform(M, P, 1, 0, x).simplify_rational()
        M = balance_transform(M, P, 1, 0, x).simplify_rational()
        M = balance_transform(M, P, 1, 0, x).simplify_rational()
        M = balance_transform(M, P, 1, 0, x).simplify_rational()

        MM, T = fuchsify(M, x)
        MM = MM.simplify_rational()
        t.assertEqual(MM, transform(M, x, T).simplify_rational())

        pranks = singularities(MM, x).values()
        t.assertEqual(pranks, [0]*len(pranks))

    def test_factorize_1(t):
        x = SR.var("x")
        e = SR.var("epsilon")
        M = matrix([[1/x, 0, 0], [0, 2/x, 0], [0, 0, 3/x]])*e
        M = transform(M, x, matrix([[1, 1, 0], [0, 1, 0],[1+2*e, 0, e]]))
        F,T = factorize(M, x, e)
        F = F.simplify_rational()
        for f in F.list():
            t.assertEqual(limit_fixed(f, e, 0), 0)

    def test_factorize_2(t):
        x = SR.var("x")
        e = SR.var("epsilon")
        M = matrix([[e*e/x]])
        with t.assertRaises(FuchsiaError):
            F,T = factorize(M, x, e)

    def test_simplify_by_jordanification(t):
        x = SR.var("x")
        M = matrix([
            [4/(x+1), -1/(6*x*(x+1)), -1/(3*x*(x+1))],
            [6*(13*x+6)/(x*(x+1)),-5*(x+3)/(3*x*(x+1)),2*(x-6)/(3*x*(x+1))],
            [-63*(x-1)/(x*(x+1)), (5*x-9)/(6*x*(x+1)), -(x-18)/(3*x*(x+1))]
        ]).simplify_rational()
        MM, T = simplify_by_jordanification(M, x)
        MM = MM.simplify_rational()
        t.assertEqual(MM, transform(M, x, T).simplify_rational())
        t.assertLess(matrix_complexity(MM), matrix_complexity(M))
