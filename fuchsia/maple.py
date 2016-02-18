from sage.interfaces.maple import maple

from fuchsia.expression import new_Expression
from fuchsia.matrix import matrix

def super_reduce(m, point):
    data = ', '.join(str(ex) for ex in m.list())
    ms = "Matrix(%d, %d, [%s])" % (m.ncols(), m.nrows(), data)
    mms = maple("convert(DEtools[super_reduce](%s, x, %s, u, 'T', 'invT')[1], list);" % (ms, point))
    data = [new_Expression(s) for s in str(mms)[1:-1].split(',')]
    mm = matrix(data, m.nrows(), m.ncols())
    t = maple.get('T')
    return mm, t
