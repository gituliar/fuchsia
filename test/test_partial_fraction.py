import unittest
import random

from   sage.all import SR, I
import fuchsia
import time

class Test(unittest.TestCase):
    def assert_pf(t, ex, x):
        pf = SR(0)
        for pi, ki, ci in fuchsia.partialer_fraction(ex, x):
            if ki >= 0:
                t.assertEqual(pi, 0)
            pf += ci*(x-pi)**ki
        t.assertEqual(ex.simplify_rational(), pf.simplify_rational())

    def test_partialer_fraction_0(t):
        x = SR.var("x")
        ex = SR(1)
        t.assert_pf(ex, x)

    def test_partialer_fraction_1(t):
        x = SR.var("x")
        ex = x**2 + x + 1 + x**-1 + x**-2
        t.assert_pf(ex, x)

    def test_partialer_fraction_2(t):
        x = SR.var("x")
        ex = 1/x/(x-1)/(x-2)
        t.assert_pf(ex, x)

    def test_partialer_fraction_3(t):
        x = SR.var("x")
        ex = 1/((x-1/(x-3)**2)-x+2)/(x-1)**2
        t.assert_pf(ex, x)

    def test_partialer_fraction_4(t):
        x = SR.var("x")
        ex = 13/((x+1)*(x-1)**3*x**3)
        t.assert_pf(ex, x)

    def test_partialer_fraction_5(t):
        x, eps = SR.var("x eps")
        ex = \
                -36*eps**3/((x+1)**3*(x-1)*x) \
                -144*eps**3/((x+1)**3*(x-1)*x**2) \
                +54*eps**2/((x+1)**3*(x-1)*x) \
                -36*eps**3/((x+1)**3*(x-1)*x**3) \
                +216*eps**2/((x+1)**3*(x-1)*x**2) \
                -26*eps/((x+1)**3*(x-1)*x) \
                +54*eps**2/((x+1)**3*(x-1)*x**3) \
                -104*eps/((x+1)**3*(x-1)*x**2) \
                +4/((x+1)**3*(x-1)*x) \
                -26*eps/((x+1)**3*(x-1)*x**3) \
                +16/((x+1)**3*(x-1)*x**2) \
                +4/((x+1)**3*(x-1)*x**3)
        t.assert_pf(ex, x)

    def test_partialer_fraction_6(t):
        x, eps = SR.var("x eps")
        ex = \
            +(((-18480*I)/557)*(-1671*I+52*(1671)**(1/2)))/(87*I+(1671)**(1/2)-(60*I)*x)**2 \
            +(18480*(1671+(52*I)*(1671)**(1/2))*eps)/(557*(87*I+(1671)**(1/2)-(60*I)*x)**2) \
            -(30*(-17824*I+1921*(1671)**(1/2))*eps)/(557*(87*I+(1671)**(1/2)-(60*I)*x)) \
            +(((18480*I)/557)*(1671*I+52*(1671)**(1/2)))/(-87*I+(1671)**(1/2)+(60*I)*x)**2 \
            +(18480*(1671-(52*I)*(1671)**(1/2))*eps)/(557*(-87*I+(1671)**(1/2)+(60*I)*x)**2) \
            -(30*(17824*I+1921*(1671)**(1/2))*eps)/(557*(-87*I+(1671)**(1/2)+(60*I)*x)) \
            -(49*eps)/(2*(-1+x)) \
            +(101*eps)/x \
            +(57*eps)/(-7+2*x)
        t.assert_pf(ex, x)

    def test_partialer_fraction_7(t):
        y, eps = SR.var("y eps")
        ex=\
                (eps - 1)*y/(y**2 + 1) - \
                2*(2*(eps - 1)*y**3 - \
                (eps - 1)*y)/(y**4 - y**2 + 1) + \
                3/2*(eps - 1)/(y + 1) + \
                3/2*(eps - 1)/(y - 1)
        t.assert_pf(ex, y)
