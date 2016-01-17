import logging

import sage.all
import sage.calculus.var
from   sage.matrix.matrix_symbolic_dense import Matrix_symbolic_dense
from   sage.misc.parser import Parser
from   sage.symbolic.expression import Expression
from   sage.symbolic.operators import mul_vararg
from   sage.symbolic.ring import SR

log = logging.getLogger('delirium')

_parser = Parser(make_var=sage.calculus.var.var)

def alphabet(ex, x):
    """Return bases of all negative exponents in the `x'-polynomial `ex'.

    >>> ex = new_Expression('a/x + b/x^2 + c/(1+x)^2')
    >>> al = alphabet(ex, 'x')
    >>> from pprint import pprint; pprint(al)
    set([x, x + 1])"""
    x = var(x)
    result = set()
    a, n = SR.wild(1), SR.wild(2)
    p_exp = a**n
    pts = singularities(ex, x)
    for pt in pts:
        r = pt.match(p_exp)
        if r:
            result.add(r[a])
        else:
            result.add(pt)
    return result

def is_Expression(obj):
    return type(obj) is Expression

def new_Expression(obj):
    if type(obj) == str:
        try:
            ex = new_Expression_from_string(obj)
        except SyntaxError:
            log.error(obj)
            raise
    else:
        raise NotImplementedError
    return ex

def new_Expression_from_string(s):
    ex = _parser.parse(s)
    return ex

def singularities(ex, x):
    x = var(x)
    result = set()
    d = ex.denominator()
    ds = [d] if d.operator() is not mul_vararg else d.operands()
    for d in ds:
        if x not in d.variables():
            continue
        result.add(d)
    return result

def var(x):
    return sage.calculus.var.var(x) if type(x) == str else x
