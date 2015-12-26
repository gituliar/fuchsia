import pytest

from   delirium.expression import alphabet, new_Expression, var

x = var("x")

def test_alphabet():
    ex = new_Expression("1/x^2 + 1/(1+x)")
    assert alphabet(ex, x) == set([x, x+1])
