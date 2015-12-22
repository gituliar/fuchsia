from   sage.calculus.var import var
from   sage.matrix.constructor import matrix
from   sage.matrix.matrix_symbolic_dense import Matrix_symbolic_dense
from   sage.misc.parser import Parser
from   sage.symbolic.expression import Expression
from   sage.symbolic.operators import mul_vararg
from   sage.symbolic.ring import SR

class Notebook(object):

    _parser = Parser(make_var=var)

    def alphabet(self, obj, x):
        """Determine the alphabet of the object, i.e. an expression or matrix.
        An alphabet is a set of x-dependent polynomials raised to a negative
        integer power.

        >>> nb = Notebook()
        >>> ex = nb.new_Expression_from_string('a/x + b/x^2 + c/(1+x)^2')
        >>> nb.alphabet(ex, 'x')
        set([x, x + 1])"""
        if type(x) == str:
            x = var(x)
        result = set()
        if self.is_Expression(obj):
            a, n = SR.wild(1), SR.wild(2)
            p_exp = a**n
            pts = self.singularities(obj, x)
            for pt in pts:
                r = pt.match(p_exp)
                if r:
                    result.add(r[a])
                else:
                    result.add(pt)
        elif self.is_Matrix(obj):
            [result.update(self.alphabet(expr,x)) for expr in obj.list()]
        return result

    def is_Expression(self, obj):
        return type(obj) is Expression

    def is_Matrix(self, obj):
        return type(obj) is Matrix_symbolic_dense

    def matrix_coefficient(self, m, x):
        """Return a leading coefficient of the matrix `m' at the point `x'."""
        x = var(x) if type(x) == str else x
        r = self.rank(m, x)
        data = []
        for expr in m.list():
            #expr = expr.partial_fraction(x)
            c = expr.coefficient(x, r).substitute(x=0)
            data.append(c)
        result = self.new_Matrix_from_list(m.ncols(), m.nrows(), data)
        return (r, result)

    def new_Expression_from_string(self, string):
        expr = self._parser.parse(string)
        return expr

    def new_Matrix_from_file(self, f):
        ncol, nrow = map(int, f.readline().split())
        try:
            data = [self.new_Expression_from_string(s) for s in f.readlines()]
        except SyntaxError:
            raise
        m = matrix(ncol, nrow, data)
        return m

    def new_Matrix_from_list(self, ncol, nrow, data):
        return matrix(ncol, nrow, data)

    def rank(self, m, x):
        """Return a Poincare rank of the matrix `m' at the point `x'."""
        x = var(x) if type(x) == str else x
        vs = [e.low_degree(x) for e in m.list() if self.is_Expression(e)]
        result = min(vs)
        return result

    def singularities(self, obj, x):
        if type(x) is str:
            x = var(x)
        result = set()
        if self.is_Expression(obj):
            d = obj.denominator()
            ds = [d] if d.operator() is not mul_vararg else d.operands()
            for d in ds:
                if x not in d.variables():
                    continue
                result.add(d)
        elif self.is_Matrix(obj):
            [result.update(self.singularities(expr,x)) for expr in obj.list()]
        return result
