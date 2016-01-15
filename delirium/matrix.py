import logging

import sage.all
from   sage.matrix.constructor import matrix as sage_matrix

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
    res = m.apply_map(lambda ex: ex.coefficient(x, n).substitute(x=0))
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

def matrix(data, nrows=None, ncols=None):
    if (nrows and ncols) is None:
        m = sage_matrix(data)
    else:
        m = sage_matrix(ncols, nrows, data)
    m = m.transpose()
    return m

def export_matrix(f, m):
    f.write("%%MatrixMarket matrix array symbolic general\n")
    f.write("%d %d\n%%\n" % (m.nrows(), m.ncols()))
    for col in m.columns():
        for mij in col:
            f.write(str(mij))
            f.write('\n')
        f.write('%\n')

def import_matrix(f):
    """Read a matrix from the file in the Matrix Market array format."""
    while True:
        s = f.readline()
        if not s.startswith('%'):
            break
    nrows, ncols = map(int, s.split())

    try:
        data = [new_Expression(s) for s in f.readlines() if not s.startswith('%')]
    except SyntaxError:
        raise

    m = matrix(data, nrows, ncols)
    return m

def fuchsian_form(m, x):
    if has_fuchsian_form_at_point(m, x):
        log.info("\nThe matrix already has a Fuchsian form.")
        return None, None
    x, n = var(x), dims(m)[0]
    m0 = coefficient_low(m, x)
    log.info("\nIs A_0 nilpotent? " + str(is_nilpotent(m0)))
    log.info("\nA_0 =\n" + m0.str(unicode=True))
    mj, tj = jordan_form(m0, x)
    mm, tm = m, matrix([0]*n**2)
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

def partial_fraction(m, x):
    res = m.apply_map(lambda obj: obj.partial_fraction(x) if is_Expression(obj) else obj)
    return res

def zero_cols(m):
    res, i = [], 1
    for col in m.columns():
        if is_zero(col):
            res.append(i)
        i += 1

def zero_rows(m):
    pass
