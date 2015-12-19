import pytest

from   sage.calculus.var import var

from   delirium.io import matrix_read
from   delirium.matrix import alphabet

def test_alphabet():
    eps, x = var("eps, x")
    with open('test/data/git_410') as f:
        m = matrix_read(f)
    assert alphabet(m, x) == set([x-1, x, x+1])

