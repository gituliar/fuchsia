import pytest

import sage.all
from   sage.calculus.var import var
from   sage.matrix.constructor import matrix
from   sage.symbolic.ring import SR

from   delirium.notebook import Notebook

nb = Notebook()
eps, x = var("eps, x")

def test_alphabet():
    with open('test/data/git_410') as f:
        m = nb.new_Matrix_from_file(f)
    assert nb.alphabet(m, x) == set([x-1, x, x+1])

def test_matrix_coefficient():
    with open('test/data/henn_324') as f:
        m = nb.new_Matrix_from_file(f)
    assert nb.matrix_coefficient(m, 'x') == (-2, matrix(SR, [[0, 0], [-1, 0]]))

def test_valuation():
    with open('test/data/henn_324') as f:
        m = nb.new_Matrix_from_file(f)
    assert nb.valuation(m, 'x') == -2
    with open('test/data/bar_ex1') as f:
        m = nb.new_Matrix_from_file(f)
    assert nb.valuation(m, 'x') == -9

def test_new_Matrix_from_file():
    with open('test/data/henn_324') as f:
        m = nb.new_Matrix_from_file(f)
    assert m == matrix(SR, [[eps/x, 0], [-1/x**2, eps/(1+x)]])
