import unittest
import random

from   sage.all import SR, matrix, oo, sqrt, I
import fuchsia

class Test(unittest.TestCase):
    def test_matrix_1(t):
        x = SR.var("x")
        M = matrix([[x**2+2+3/(x-1)**2]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        t.assertTrue((m1.get_M() - M).is_zero())
        t.assertTrue((m2.get_M() - M).is_zero())

    def test_matrix_2(t):
        x = SR.var("x")
        M = matrix([[1/x+2/(x-1)**2+3+4*x+5*x*x, 1/(x-1)/(x-2)], [x*(x-1), (x+1)/(x+2)]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        t.assertTrue((m1.get_M() - M).is_zero())
        t.assertTrue((m2.get_M() - M).is_zero())

    def test_matrix_3(t):
        return
        x, eps = SR.var("x eps")
        M = matrix([
            [eps/(x - 1), 0, 0],
            [(x**3 + x**2 - 3*x - 3)/(2*x**3 - 2*x**2 - x + 1), 2*eps/(x - 3), 0],
            [-(x**3 + 3*x**2 - 3*x + 2)/(2*x**3 + x**2 - 2*x) + 2*(x**3 - x**2 - x + 1)/(2*x + 1),
                -3*(x**3 + x - 1)/(3*x**3 + x - 2),
                4*eps/(x - 5)]
        ])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        t.assertTrue((m1.get_M() - M).is_zero())
        t.assertTrue((m2.get_M() - M).is_zero())

    def test_singularities_1(t):
        x = SR.var("x")
        M = matrix([[1/x+2/(x-1)**2+3+4*x+5*x*x, 1/(x-1)/(x-2)], [x*(x-1), (x+1)/(x+2)]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        t.assertTrue((m1.get_M() - M).is_zero())
        t.assertTrue((m2.get_M() - M).is_zero())
        t.assertEqual(m1.singular_points(), m2.singular_points())

    def test_singularities_2(t):
        x = SR.var("x")
        M = matrix([[1/x+1/(x-1)]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        t.assertTrue((m1.get_M() - M).is_zero())
        t.assertTrue((m2.get_M() - M).is_zero())
        t.assertEqual(m1.singular_points(), m2.singular_points())

    def test_singularities_3(t):
        x = SR.var("x")
        M = matrix([[1/x-1/(x-1)]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        t.assertTrue((m1.get_M() - M).is_zero())
        t.assertTrue((m2.get_M() - M).is_zero())
        t.assertEqual(m1.singular_points(), m2.singular_points())

    def test_c0_and_c1_1(t):
        x = SR.var("x")
        M = matrix([[1/(x-1)+5/(x-2)+7/(x-3)]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        for p, prank in m1.singular_points().iteritems():
            t.assertEqual(m1.c0(p, prank), m2.c0(p, prank))
            #t.assertTrue((m1.c0(p, prank) - m2.c0(p, prank)).is_zero())
            if prank > 0:
                t.assertEqual(m1.c1(p, prank), m2.c1(p, prank))
                #t.assertTrue((m1.c1(p, prank) - m2.c1(p, prank)).is_zero())

    def test_c0_and_c1_2(t):
        x = SR.var("x")
        M = matrix([[x+2+3/(x-1)+4/(x-2)**2]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        for p, prank in m1.singular_points().iteritems():
            t.assertEqual(m1.c0(p, prank), m2.c0(p, prank))
            #t.assertTrue((m1.c0(p, prank) - m2.c0(p, prank)).is_zero())
            if prank > 0:
                t.assertEqual(m1.c1(p, prank), m2.c1(p, prank))
                #t.assertTrue((m1.c1(p, prank) - m2.c1(p, prank)).is_zero())

    def test_c0_and_c1_3(t):
        x = SR.var("x")
        M = matrix([[2+3/(x-1)**2+4/(x-2)**2, 1/x-2/x**2], [x+2*x*x, 1/(x*x-2)]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        for p, prank in m1.singular_points().iteritems():
            t.assertEqual(m1.c0(p, prank), m2.c0(p, prank))
            #t.assertTrue((m1.c0(p, prank) - m2.c0(p, prank)).is_zero())
            if prank > 0:
                t.assertEqual(m1.c1(p, prank), m2.c1(p, prank))
                #t.assertTrue((m1.c1(p, prank) - m2.c1(p, prank)).is_zero())

    def test_balance_1(t):
        x = SR.var("x")
        M = matrix([[1/(x-2)+4/(x-5)]])
        #M = matrix([[11/(x-5)+7/(x-2)**4+5/(x-2)**3+3/(x-2)**2+1/(x-2)]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        P = matrix([[1]], ring=SR)
        t.assertTrue((P*P-P).is_zero())
        def testB(x1, x2):
            m1b = m1.apply_balance(P, x1, x2)
            m2b = m2.apply_balance(P, x1, x2)
            t.assertEqual(m1b.get_T(), m2b.get_T())
            t.assertEqual(
                    fuchsia.partial_fraction(m1b.get_M(), x),
                    fuchsia.partial_fraction(m2b.get_M(), x))
        testB(2, 5)
        testB(5, 2)
        testB(2, 1)
        testB(1, 5)
        testB(-3, 7)

    def test_balance_2(t):
        x = SR.var("x")
        M = matrix([
            [10/x, 11/x],
            [0, 12/x]
        ])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        P = matrix([[1,0],[1,0]], ring=SR)
        t.assertTrue((P*P-P).is_zero())
        def testB(x1, x2):
            m1b = m1.apply_balance(P, x1, x2)
            m2b = m2.apply_balance(P, x1, x2)
            t.assertEqual(m1b.get_T(), m2b.get_T())
            t.assertEqual(
                    fuchsia.partial_fraction(m1b.get_M(), x),
                    fuchsia.partial_fraction(m2b.get_M(), x))
        testB(0, 5)
        testB(5, 0)
        testB(2, 1)
        testB(-3, 7)

    def test_balance_3(t):
        x = SR.var("x")
        M = matrix([[1/x/x, 2/(x-1)], [3/x, 4/(x-1)**2]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        P = matrix([[1, 0], [0, 1]], ring=SR)
        def testB(x1, x2):
            m1b = m1.apply_balance(P, x1, x2)
            m2b = m2.apply_balance(P, x1, x2)
            t.assertEqual(m1.get_T(), m2.get_T())
            t.assertEqual(
                    fuchsia.partial_fraction(m1b.get_M(), x),
                    fuchsia.partial_fraction(m2b.get_M(), x))
        testB(0, 1)
        testB(1, 0)
        testB(0, 3)
        testB(3, 1)
        testB(-3, 7)

    def test_balance_4(t):
        x = SR.var("x")
        M = matrix([[x+2, 3/(x-1)], [3*x, 4/(x-1)**2]])
        m1 = fuchsia.GeneralSystem.from_M(M, x)
        m2 = fuchsia.RationalSystem.from_M(M, x)
        P = matrix([[1, 0], [0, 1]], ring=SR)
        def testB(x1, x2):
            m1b = m1.apply_balance(P, x1, x2)
            m2b = m2.apply_balance(P, x1, x2)
            t.assertEqual(m1.get_T(), m2.get_T())
            t.assertEqual(
                    fuchsia.partial_fraction(m1b.get_M(), x),
                    fuchsia.partial_fraction(m2b.get_M(), x))
        testB(oo, 1)
        testB(1, oo)
        testB(1, oo)
        testB(oo, 3)
        testB(3, 1)
        testB(-3, 7)

    def test_partialer_fraction_1(t):
        x = SR.var("x")
        t.assertEqual(fuchsia.partialer_fraction(SR(0), x), [])

    def test_partialer_fraction_2(t):
        x = SR.var("x")
        t.assertEqual(
                fuchsia.partialer_fraction(SR(9), x),
                [(SR(0), 0, SR(9))])

    def test_partialer_fraction_3(t):
        x = SR.var("x")
        ex = 3*x**2 + 4*x + 5 + 6/(x - 2)**3
        t.assertEqual(
                set(fuchsia.partialer_fraction(ex, x)),
                set([
                    (SR(0), 0, SR(5)),
                    (SR(0), 1, SR(4)),
                    (SR(0), 2, SR(3)),
                    (SR(2), -3, SR(6))
                ]))

    def test_partialer_fraction_4(t):
        x = SR.var("x")
        ex = (x*3+1)/(x**3+x**2+x+1)+1/(x**2-5)+(x**3*(2+3/x)+2)/(x**2+7)
        pf = SR(0)
        for p, k, c in fuchsia.partialer_fraction(ex, x):
            t.assertEqual(p.variables(), ())
            t.assertTrue(type(k) is int)
            if k >= 0: t.assertEqual(p, SR(0))
            t.assertFalse(c.is_zero())
            t.assertEqual(c.variables(), ())
            pf += c*(x-p)**k
        t.assertEqual(pf, ex)
