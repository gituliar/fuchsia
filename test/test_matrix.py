import pytest

from   delirium.matrix import (
           alphabet, coefficient_low, degree_high, degree_low, new_Matrix, var)

eps, x = var("eps, x")

def test_alphabet():
    with open('test/data/git_410') as f:
        m = new_Matrix(f)
    assert alphabet(m, x) == set([x-1, x, x+1])

def test_coefficient_low():
    with open('test/data/henn_324') as f:
        m = new_Matrix(f)
    assert coefficient_low(m, 'x') == new_Matrix([0, 0, -1, 0])

def test_degree():
    # test 01
    with open('test/data/henn_324') as f:
        m = new_Matrix(f)
    assert degree_high(m, 'x') == 0
    assert degree_low(m, 'x') == -2
    # test 02
    with open('test/data/bar_ex1') as f:
        m = new_Matrix(f)
    assert degree_high(m, 'x') == 0
    assert degree_low(m, 'x') == -9

def test_new_Matrix_from_file():
    with open('test/data/henn_324') as f:
        m = new_Matrix(f)
    assert m == new_Matrix([eps/x, 0, -1/x**2, eps/(1+x)])
