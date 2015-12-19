from   sage.calculus.var import var
from   sage.matrix.constructor import matrix
from   sage.misc.parser import Parser
from   sage.symbolic.ring import SR

def matrix_read(f):
    """Read a matrix from the open file"""
    nc, nr = map(int, f.readline().split())
    m = matrix(SR, nr, nc)
    p = Parser(make_var=var)
    for i in xrange(nr):
        for j in xrange(nr):
            s = f.readline()
            try:
                e = p.parse(s)
            except SyntaxError:
                print s
                return None
            m[i,j] = e
    return m
