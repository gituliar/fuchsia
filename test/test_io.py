import pytest

import sage.all
from   sage.calculus.var import var
from   sage.matrix.constructor import matrix
from   sage.symbolic.ring import SR

from   delirium.io import matrix_read

def test_tokenizerer_wolfram():
    eps, x = var("eps, x")
    with open('test/data/henn_324') as f:
        m = matrix_read(f)
    assert m == matrix(SR, [[eps/x, 0], [-1/x**2, eps/(1+x)]])
