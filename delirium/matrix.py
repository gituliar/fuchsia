import logging

from   sage.matrix.constructor import matrix

from   delirium.expression import is_Expression, new_Expression, var
from   delirium.expression import alphabet as ex_alphabet

log = logging.getLogger('delirium')

def _matrix_list_map(m, f):
    res = [f(item) for item in m.list()]
    return res


def alphabet(m, x):
    res = set()
    _matrix_list_map(m, lambda ex: res.update(ex_alphabet(ex,x)))
    return res

def coefficient(m, x, n):
    """Return a coefficient of the exponent x^n in the matrix `m'."""
    x = var(x)
    data = _matrix_list_map(m, lambda ex: ex.coefficient(x, n).substitute(x=0))
    res = new_Matrix(data)
    return res

def coefficient_low(m, x, n=0):
    d = n + degree_low(m, x)
    res = coefficient(m, x, d)
    return res

def degree_high(m, x):
    """Return the highest integer exponent of the base `x' in the matrix `m'."""
    x = var(x)
    ds = _matrix_list_map(m, lambda ex: ex.degree(x) if is_Expression(ex) else None)
    res = max(ds)
    return res

def degree_low(m, x):
    """Return the lowest integer exponent of the base `x' in the matrix `m'."""
    x = var(x)
    ds = _matrix_list_map(m, lambda ex: ex.low_degree(x) if is_Expression(ex) else None)
    res = min(ds)
    return res

def dims(m):
    """Return matrix size as (cols, rows)."""
    return m.dimensions()

def fuchsian_form(m, x):
    if has_fuchsian_form_at_point(m, x):
        log.info("\nThe matrix already has a Fuchsian form.")
        return None, None
    x, n = var(x), dims(m)[0]
    m0 = coefficient_low(m, x)
    log.info("\nIs A_0 nilpotent? " + str(is_nilpotent(m0)))
    log.info("\nA_0 =\n" + m0.str(unicode=True))
    mj, tj = jordan_form(m0, x)
    mm, tm = m, new_Matrix([0]*n**2)
    return mm, tm

def has_fuchsian_form(m, x):
    res = True
    for letter in alphabet(m, x):
        if not has_fuchsian_at_point(m, letter):
            res = False
            break
    return res

def has_fuchsian_form_at_point(m, x):
    return degree_low(m, x) >= -1

def is_Matrix(obj):
    return type(obj) is Matrix_symbolic_dense

def is_nilpotent(m):
    res = is_zero(m.eigenvalues())
    return res

def is_zero(obj):
    if type(obj) in (list, tuple):
        res = True
        for item in obj:
            if item != 0:
                res = False
                break
    elif is_Matrix(obj):
        res = is_zero(obj.list())
    else:
        raise NotImplementedError
    return res

def jordan_form(m, x):
    x = var(x)
    mj, tj = m.jordan_form(transformation=True)
    log.info("\nJordan form =\n" + str(mj))
    log.info("\nJordan transformation =\n" + str(tj))
    return mj, tj

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

def zero_cols(m):
    res, i = [], 1
    for col in m.columns():
        if is_zero(col):
            res.append(i)
        i += 1

def zero_rows(m):
    pass
