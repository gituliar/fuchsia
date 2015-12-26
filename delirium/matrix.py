from   sage.matrix.constructor import matrix

from   delirium.expression import is_Expression, new_Expression, var
from   delirium.expression import alphabet as ex_alphabet


def _matrix_list_map(m, x, f):
    x = var(x)
    res = [f(item) for item in m.list()]
    return res


def alphabet(m, x):
    res = set()
    _matrix_list_map(m, x, lambda ex: res.update(ex_alphabet(ex,x)))
    return res

def coefficient(m, x, n=None):
    """Return the coefficient in `x' the matrix `m' at the point `x' and the valuation `n'.
    In other words, return the matrix m_n, such that m = sum_n m_n x^n."""
    x = var(x)
    d = n if n is not None else degree_low(m, x)
    data = _matrix_list_map(m, x, lambda ex: ex.coefficient(x, d).substitute(x=0))
    res = new_Matrix(data)
    return res

def degree_high(m, x):
    """Return the highest integer exponent of the base `x' in the matrix `m'."""
    x = var(x)
    ds = _matrix_list_map(m, x, lambda ex: ex.degree(x) if is_Expression(ex) else None)
    res = max(ds)
    return res

def degree_low(m, x):
    """Return the lowest integer exponent of the base `x' in the matrix `m'."""
    x = var(x)
    ds = _matrix_list_map(m, x, lambda ex: ex.low_degree(x) if is_Expression(ex) else None)
    res = min(ds)
    return res

def is_Matrix(obj):
    return type(obj) is Matrix_symbolic_dense

def new_Matrix(obj):
    if type(obj) is list:
        res = new_Matrix_from_list(obj)
    elif type(obj) is file:
        res = new_Matrix_from_file(obj)
    else:
        raise NotImplementedError
    return res

def new_Matrix_from_file(f):
    ncol, nrow = map(int, f.readline().split())
    try:
        data = [new_Expression(s) for s in f.readlines()]
    except SyntaxError:
        raise
    m = new_Matrix_from_list(data)
    return m

def new_Matrix_from_list(data):
    n = len(data)
    ncol = int(n**0.5)
    assert n == ncol**2
    m = matrix(ncol, ncol, data)
    return m
