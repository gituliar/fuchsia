#!/usr/bin/env sage
"""\
Usage:
    fuchsia [options] <command> <args>...

Commands:
    reduce [-x <name>] [-e <name>] [-m <path>] [-t <path>] <matrix>
        find an epsilon form of the given matrix

    fuchsify [-x <name>] [-m <path>] [-t <path>] <matrix>
        find a transformation that will transform a given matrix
        into Fuchsian form

    normalize [-x <name>] [-e <name>] [-m <path>] [-t <path>] <matrix>
        find a transformation that will transform a given Fuchsian
        matrix into normalized form

    factorize [-x <name>] [-e <name>] [-m <path>] [-t <path>] <matrix>
        find a transformation that will make a given normalized
        matrix proportional to the infinitesimal parameter

    info [-x <name>] <matrix>
        show a description of a given matrix, listing singular
        points and residue eigenvalues (in verbose more)

    sort [-m <path>] [-t <path>] <matrix>
        find a block-triangular form of the given matrix

    transform [-x <name>] [-m <path>] <matrix> <transform>
        transform a given matrix using a given transformation

    changevar [-x <name>] [-y <name>] [-m <path>] <matrix> <expr>
        transform a given matrix by susbtituting free variable
        by a given expression

Options:
    -h          show this help message
    -f <fmt>    matrix file format: mtx or m (default: m)
    -l <path>   write log to this file
    -v          produce a more verbose log
    -P <path>   save profile report into this file
    -x <name>   use this name for the free variable (default: x)
    -y <name>   use this name for the new free variable (default: y)
    -e <name>   use this name for the infinitesimal parameter (default: eps)
    -m <path>   save the resulting matrix into this file
    -t <path>   save the resulting transformation into this file
    -s <number> use this random seed when applicable (default: 0)
    --use-maple speed up calculations by using Maple when possible

Arguments:
    <matrix>    read the input matrix from this file
    <transform> read the transformation matrix from this file
    <expr>      arbitrary expression
"""

__author__ = "Oleksandr Gituliar, Vitaly Magerya"
__author_email__ = "oleksandr@gituliar.net"
__version__ = "17.06.28"

__all__ = [
    "FuchsianSystem",
    "NormalizeAssistant",

    "balance",
    "balance_transform",
    "block_triangular_form",
    "epsilon_form",
    "export_matrix_to_file",
    "factorize",
    "fuchsify",
    "fuchsify_off_diagonal_blocks",
    "import_matrix_from_file",
    "is_fuchsian",
    "matrix_c0",
    "matrix_c1",
    "matrix_complexity",
    "matrix_residue",
    "normalize",
    "reduce_diagonal_blocks",
    "setup_fuchsia",
    "simple",
    "simplify_by_factorization",
    "simplify_by_jordanification",
    "singularities",
    "transform"
]

from   collections import defaultdict
from   functools import wraps
from   itertools import combinations, permutations
from   random import Random
import inspect
import logging
import time

from   sage.all import *
from   sage.misc.parser import Parser
from   sage.libs.ecl import ecl_eval

def logcall(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        logger.enter(f.__name__)
        r = f(*args, **kwargs)
        logger.exit(f.__name__)
        return r
    return wrapper

class ElapsedTimeFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt=None):
        super(ElapsedTimeFormatter, self).__init__(fmt, datefmt)
        self.start_time = time.time()
    def formatTime(self, record, datefmt=None):
        if datefmt is None:
            return "%7.1fs" % (record.created - self.start_time)
        return time.strftime(datefmt, time.localtime(record.created))

class FuchsiaLogger(object):
    def __init__(self):
        log_handler = logging.StreamHandler()
        log_handler.setFormatter(ElapsedTimeFormatter(
            "\033[32m[%(asctime)s]\033[0m %(message)s"
        ))
        logger = logging.getLogger('fuchsia')
        logger.addHandler(log_handler)
        logger.setLevel(logging.WARNING)

        self.addHandler = logger.addHandler
        self.error = logger.error
        self.isEnabledFor = logger.isEnabledFor
        self.logger = logger
        self.setLevel = logger.setLevel

        self.depth = 0

    def indent(self):
        return "   " * self.depth

    def is_verbose(self):
        return self.isEnabledFor(logging.INFO)

    def enter(self, name):
        self.info("-> %s" % name)
        self.depth += 1

    def exit(self, name):
        self.depth -= 1
        self.info("<- %s" % name)

    def debug(self, msg):
        self.logger.debug(self.indent() + msg)

    def info(self, msg):
        self.logger.info(self.indent() + msg)

if True:
    ecl_eval("(ext:set-limit 'ext:heap-size 0)")
    logger = FuchsiaLogger()
    USE_MAPLE = False

def setup_fuchsia(verbosity=0, use_maple=False):
    global USE_MAPLE
    USE_MAPLE = bool(use_maple)
    logger.setLevel(
        logging.WARNING if verbosity <= 0 else \
        logging.INFO if verbosity == 1 else
        logging.DEBUG
    )

class FuchsiaError(Exception):
    pass

def cross_product(v1, v2):
    m1, m2 = matrix(v1), matrix(v2)
    return m1.transpose() * m2

def dot_product(v1, v2):
    m1, m2 = matrix(v1), matrix(v2)
    sp = m1 * m2.transpose()
    return sp[0,0]

def partial_fraction(M, var):
    return M.apply_map(lambda ex: ex.partial_fraction(var))

#@logcall
def fuchsia_simplify(obj, x=None):
    if USE_MAPLE:
        def maple_simplify(ex):
            if hasattr(ex, "variables") and ((x is None) or (x in ex.variables())):
                res = maple.factor(maple.radnormal(ex))
                return parse(str(res))
            else:
                return ex
        if hasattr(obj, "apply_map"):
            return obj.apply_map(maple_simplify)
        else:
            return maple_simplify(obj)
    else:
        return obj.simplify_rational()

def fuchsia_solve(eqs, var):
    if USE_MAPLE:
        s = maple.solve(eqs, var)
        solutions = s.parent().get(s._name).strip('[]').split('],[')
        solutions = [s.split(',') for s in solutions if s != '']
        result = []
        for solution in solutions:
            r = []
            for s in solution:
                try:
                    expr = fuchsia_simplify(parse(s))
                    r.append(expr)
                except SyntaxError as error:
                    print("ERROR:  \n%s\n  %s\n" % (s, error))
                    continue
            result.append(r)
        return result
    else:
        return solve(eqs, var, solution_dict=True)

def change_variable(m, x, y, fy):
    mm = m.subs({x: fy}) * derivative(fy, y)
    return mm

def transform(M, x, T):
    """Given a system of differential equations dF/dx = M*F,
    and a transformation of base functions F = T*F', compute
    and return M', such that dF'/dx = M'*F'.

    Note: M' = inverse(T)*(M*T - dT/dx)
    """
    mm = T.inverse()*(M*T - derivative(T, x))
    return mm

def balance(P, x1, x2, x):
    assert P.is_square()
    assert (P*P - P).is_zero()
    assert not bool(x1 == x2)
    coP = identity_matrix(P.nrows()) - P
    if x1 == oo:
        b = coP - (x - x2)*P
    elif x2 == oo:
        b = coP - 1/(x - x1)*P
    else:
        b = coP + (x - x2)/(x - x1)*P
    return b

def balance_transform(M, P, x1, x2, x):
    """Same thing as transform(M, x, balance(P, x1, x2, x)), but faster."""
    assert P.is_square()
    #assert (P*P - P).is_zero()
    assert not bool(x1 == x2)
    coP = identity_matrix(P.nrows()) - P
    if x1 == oo:
        k = -(x - x2)
        d = -1
    elif x2 == oo:
        k = -1/(x - x1)
        d = 1/(x - x1)**2
    else:
        k = (x - x2)/(x - x1)
        d = (x2 - x1)/(x - x1)**2
    mm = (coP + 1/k*P)*M*(coP + k*P) - d/k*P
    return mm

#@logcall
def limit_fixed(expr, x, x0):
    #logger.debug("[%s -> %s] %s" % (x, x0, expr))
    res = None
    if USE_MAPLE:
        res = limit_fixed_maple(expr, x, x0)
    else:
        res = limit_fixed_maxima(expr, x, x0)
    #logger.debug("res = %s" % res)
    return res

def limit_fixed_maple(expr, x, x0):
    res = maple.limit(expr, **{str(x): x0})
    try:
        res = parse(str(res))
    except:
        res = NaN
    return res

def limit_fixed_maxima(expr, x, x0):
    """Return a limit of expr when x->lim.

    The standard 'limit()' function of SageMath does not allow
    you to specify the variable, only it's name as a keyword
    argument. If you have a variable, and not it's name, use
    this function instead.
    """
    l = maxima_calculus.sr_limit(expr, x, x0)
    return expr.parent()(l)

@logcall
def singularities(m, x):
    """Find values of x around which rational matrix M has
    a singularity; return a dictionary with {val: p} entries,
    where p is the Poincare rank of M at x=val.

    Example:
    >>> x = var("x")
    >>> s = singularities(matrix([[1/x, 0], [1/(x+1)**3, 1/(x+1)]]), x)
    >>> sorted(s.items())
    [(-1, 2), (0, 0), (+Infinity, 0)]

    >>> s = singularities(matrix([[x, 1, 1/x]]), x)
    >>> sorted(s.items())
    [(0, 0), (+Infinity, 2)]
    """
    m = fuchsia_simplify(m)
    result = {}
    for expr in m.list():
        points_expr = singularities_expr(expr, x)
        for x0, p in points_expr.items():
            if x0 in result:
                result[x0] = max(result[x0], p)
            else:
                result[x0] = p
    return result

def singularities_expr(expr, x):
    if USE_MAPLE:
        return singularities_expr_maple(expr, x)
    else:
        return singularities_expr_maxima(expr, x)

def singularities_expr_maple(expr, x):
    if bool(expr == 0):
        return {}

    result = {}
    sols = fuchsia_solve(1/expr, x)
    points = [x0 for x0 in sols[0]] if len(sols) > 0 else []
    for x0 in points:
        if x0 not in result:
            result[x0] = 0
        else:
            result[x0] += 1

    sols = fuchsia_solve((1/(expr.subs({x: 1/x})/x**2)).simplify_rational(), x)
    points = [x0 for x0 in sols[0]] if len(sols) > 0 else []
    for x0 in points:
        if x0 == 0:
            if oo not in result:
                result[oo] = 0
            else:
                result[oo] += 1
    return result

def singularities_expr_maxima(expr, x):
    if expr.is_zero():
        return {}

    result = {}
    eq = factor(1/expr)
    points, ps = solve(eq, x, solution_dict=True, multiplicities=True)
    for sol, p in zip(points, ps):
        if x not in sol:
            raise Exception("Maxima can't solve this equation: %s" % eq)
        x0 = expand(sol[x])
        if p >= 1:
            result[x0] = p-1

    eq = factor(1/(expr.subs({x: 1/x})/x**2))
    points, ps = solve(eq, x, solution_dict=True, multiplicities=True)
    for sol, p in zip(points, ps):
        if x not in sol:
            raise Exception("Maxima can't solve this equation: %s" % eq)
        if sol[x] == 0 and p >= 1:
            result[oo] = p-1
    return result

def matrix_taylor0(m, x, x0, exp):
    """Return the 0-th coefficient of Taylor expansion of
    a matrix M around a finite point x=point, assuming that
    M(x->point)~1/(x-point)**exp.

    Example:
    >>> x = var('x')
    >>> matrix_taylor0(matrix([[x/(x-1), 1/x, x, 1]]), x, 1, 1)
    [1 0 0 0]
    """
    return m.apply_map(lambda ex: limit_fixed(ex*(x-x0)**exp, x, x0))

def matrix_taylor1(M, x, point, exp):
    """Return the 1-th coefficient of Taylor expansion of
    a matrix M around a finite point x=point, assuming that
    M(x->point)~1/(x-point)**exp.

    Example:
    >>> x = var('x')
    >>> matrix_taylor1(matrix([[x/(x-1), 1/x, x, 1]]), x, 0, 1)
    [0 0 0 1]
    """
    def taylor1(e, x):
        l0 = limit_fixed(e, x, 0)
        l1 = limit_fixed((e - l0)/x, x, 0)
        return l1
    return matrix([
        [taylor1(e, x) for e in row]
        for row in M.subs({x: x+point})*x**exp
    ])

def matrix_c0(M, x, point, p):
    """Return the 0-th coefficient of M's expansion at x=point,
    assuming Poincare rank of M at that point is p. If point is
    +Infinity, return minus the coefficient at the highest power of x.

    Examples:
    >>> x = var("x")
    >>> m = matrix([[1/x, 1/x**2], [1, 1/(x-1)]])
    >>> matrix_c0(m, x, 0, 1)
    [0 1]
    [0 0]
    >>> matrix_c0(m, x, 1, 0)
    [0 0]
    [0 1]
    >>> matrix_c0(m, x, oo, 1)
    [ 0  0]
    [-1  0]
    >>> matrix_c0(m*x, x, oo, 2)
    [ 0  0]
    [-1  0]
    """
    if point == oo:
        return -matrix_taylor0(M.subs({x: 1/x}), x, 0, p-1)
    else:
        return matrix_taylor0(M, x, point, p+1)

def matrix_c1(M, x, point, p):
    """Return the 1-st coefficient of M's expansion at x=point,
    assuming Poincare rank of M at that point is p. If point is
    +Infinity, return minus the coefficient at the second-to-highest
    power of x.

    Examples:
    >>> x = var("x")
    >>> m = matrix([[1/x, 1/x**2], [1, 1/(x-1)]])
    >>> matrix_c1(m, x, 0, 1)
    [1 0]
    [0 0]
    >>> matrix_c1(m, x, oo, 1)
    [-1  0]
    [ 0 -1]
    """
    if point == oo:
        return -matrix_taylor1(M.subs({x: 1/x}), x, 0, p-1)
    else:
        return matrix_taylor1(M, x, point, p+1)

def matrix_residue(M, x, x0):
    """Return matrix residue of M at x=x0, assuming that M's
    Poincare rank at x=x0 is zero.

    Example:
    >>> x = var("x")
    >>> m = matrix([[1/x, 2/x], [3/x, 4/x]])
    >>> matrix_residue(m, x, 0)
    [1 2]
    [3 4]
    >>> matrix_residue(m, x, oo)
    [-1 -2]
    [-3 -4]
    """
    if M._cache is None:
        M._cache = {}
    key = "matrix_residue_%s_%s" % (x, x0)
    if key in M._cache:
        return M._cache[key]

    m0 = matrix_c0(M, x, x0, 0)
    M._cache[key] = m0
    return m0

def matrix_is_nilpotent(M):
    """Return True if M is always nilpotent, False otherwise.

    Examples:
    >>> a, b, c = var('a b c')
    >>> matrix_is_nilpotent(matrix([[0,a,b],[0,0,c],[0,0,0]]))
    True
    >>> matrix_is_nilpotent(matrix([[a,0,0],[0,b,0],[0,0,c]]))
    False
    """
    for v in M.eigenvalues():
        if not v.is_zero():
            return False
    return True

def jordan_cell_sizes(J):
    """Return a tuple of Jordan cell sizes from a matrix J in Jordan
    normal form.

    Examples:
    >>> jordan_cell_sizes(matrix([[1,1,0,0],[0,1,0,0],[0,0,1,1],[0,0,0,1]]))
    (2, 2)
    >>> jordan_cell_sizes(matrix([[1,1,0,0],[0,1,1,0],[0,0,1,0],[0,0,0,1]]))
    (3, 1)
    >>> jordan_cell_sizes(zero_matrix(5,5))
    (1, 1, 1, 1, 1)
    """
    assert J.is_square()
    sizes = []
    n = 1
    for i in range(J.nrows() - 1):
        if J[i, i+1].is_zero():
            sizes.append(n)
            n = 1
        else:
            n += 1
    sizes.append(n)
    assert sum(sizes) == J.nrows()
    return tuple(sizes)

def solve_right_fixed(A, B):
    """As of SageMath 6.10, 'Matrix.solve_right' method uses a broken
    check for solution correctness; this function corrects that check.
    """
    C = A.solve_right(B, check=False)
    if C.nrows() and C.ncols():
        # Sometimes 'solve_right' returns a giant expression,
        # which can be simplified into a tiny one; without this
        # simplification 'is_zero' call below may take forever
        # to finish.
        #
        # Note that the condition above is there to make sure
        # that an empty matrix is not passed into 'simplify_rational',
        # otherwise you'll get this error:
        #   TypeError: unable to make sense of Maxima expression
        #   'matrix()' in Sage
        C = C.simplify_rational()
    if not (A*C - B).simplify_rational().is_zero():
        raise ValueError("matrix equation has no solutions")
    return C

def solve_left_fixed(A, B):
    """As of SageMath 6.10, 'Matrix.solve_left' method uses a broken
    check for solution correctness; this function corrects that check.
    """
    return solve_right_fixed(A.transpose(), B.transpose()).transpose()

def any_integer(rng, ring, excluded):
    r = 2
    while True:
        p = ring(rng.randint(-r, r))
        if p not in excluded:
            return p
        r *= 2

#==================================================================================================
# Transformation routines
#==================================================================================================

def block_triangular_form(m):
    """Find a lower block-triangular form of a given matrix.

    Return a tuple `(M, T, B)` where:
      * `M` is a new matrix;
      * `T` is a transformation matrix;
      * `B` is a list of tuples (ki, ni), such that M's i-th
        diagonal block is given by `M.submatrix(ki, ki, ni, ni)`.
    """
    logger.enter("block_triangular_form")
    if logger.is_verbose():
        logger.debug("matrix before transformation:\n%s\n" % matrix_mask_str(m))

    n = m.nrows()
    deps_1_1 = {}
    for i, row in enumerate(m.rows()):
        deps_1_1[i] = set(j for j,ex in enumerate(row) if not bool(ex==0))

    deps_1_all = {}
    def find_deps_1_all(i):
        if i in deps_1_all:
            return deps_1_all[i]
        deps_1_all[i] = deps_1_1[i]
        for j in deps_1_1[i]:
            if i == j:
                continue
            find_deps_1_all(j)
            deps_1_all[i] = deps_1_all[i].union(deps_1_all[j])
    [find_deps_1_all(j) for j in range(n)]

    deps_coup = dict((i, set([])) for i in range(n))
    for i in range(n):
        for j in range(n):
            if (i in deps_1_all[j]) and (j in deps_1_all[i]):
                deps_coup[i].update([i,j])
                deps_coup[j].update([i,j])

    shuffle = []
    error = False
    while not error:
        if not deps_coup:
            break
        error = True
        for i in range(n):
            if i not in deps_coup:
                continue
            b = deps_coup[i].copy()
            if b != deps_1_all[i]:
                continue
            if b == set([]):
                b = set([i])
            shuffle += [b]
            for j in deps_1_all[i].copy():
                if j in deps_1_all:
                    del deps_1_all[j]
                    del deps_coup[j]
            for j in range(n):
                if j not in deps_coup:
                    continue
                deps_coup[j].difference_update(b)
                deps_1_all[j].difference_update(b)
            if i in deps_coup:
                del deps_coup[i]
            if i in deps_1_all:
                del deps_1_all[i]
            error = False
            break
    if error:
        raise FuchsiaError("Infinite loop")

    t = zero_matrix(SR,n)
    i = 0
    blocks = []
    for block in shuffle:
        blocks.append((i, len(block)))
        for j in list(block):
            t[j,i] = 1
            i += 1
    logger.info("found %d blocks" % len(blocks))
    mt = transform(m, None, t)
    if logger.is_verbose():
        logger.debug("matrix after transformation:\n%s\n" % matrix_mask_str(mt))
    logger.exit("block_triangular_form")
    return mt, t, blocks

def epsilon_form(m, x, eps, seed=0):
    logger.enter("epsilon_form")
    m, t1, b = block_triangular_form(m)
    m, t2 = reduce_diagonal_blocks(m, x, eps, b=b, seed=seed)
    m, t3 = fuchsify_off_diagonal_blocks(m, x, eps, b=b)
    m, t4 = factorize(m, x, eps, b=b, seed=seed)
    t = t1*t2*t3*t4
    logger.exit("epsilon_form")
    return m, t

def matrix_mask(m):
    n = m.nrows()
    return matrix(SR, n, n, [int(not bool(ex==0)) for ex in m.list()])

def matrix_mask_str(m):
    s = ''
    for row in matrix_mask(m).rows():
        s += ' '.join([str(ex) for ex in row]).replace('0', '.').replace('1','x') + "\n"
    return s

def matrix_str(m, n=2):
    buf = ""
    ind = " "*n
    for col in m.columns():
        for ex in col:
            buf += ind + str(ex) + "\n"
        buf += "\n"
    buf = buf[:-1]
    return buf

#==================================================================================================
# Step I: Fuchsify
#==================================================================================================

def is_fuchsian(m, x):
    for i, expr in enumerate(m.list()):
        if expr.is_zero():
            continue
        points = singularities_expr(expr, x)
        for x0, p in points.items():
            if p > 0:
                return False
    return True

def fuchsify(M, x, seed=0):
    """Given a system of differential equations of the form dF/dx=M*F,
    try to find a transformation T, which will reduce M to Fuchsian
    form. Return the transformed M and T. Raise FuchsiaError if
    M can not be transformed into Fuchsian form.

    Note that such transformations are not unique; you can obtain
    different ones by supplying different seeds.
    """
    logger.enter("fuchsify")
    assert M.is_square()
    rng = Random(seed)
    poincare_map = singularities(M, x)
    def iter_reductions(p1, U):
        for p2, prank2 in poincare_map.items():
            if bool(p2 == p1): continue
            while prank2 >= 0:
                B0 = matrix_c0(M, x, p2, prank2)
                if not B0.is_zero(): break
                poincare_map[p2] = prank2 = prank2 - 1
            if prank2 < 0: continue
            v = find_dual_basis_spanning_left_invariant_subspace(B0, U, rng)
            if v is None: continue
            P = fuchsia_simplify(U*v, x)
            M2 = fuchsia_simplify(balance_transform(M, P, p1, p2, x), x)
            yield p2, P, M2
    combinedT = identity_matrix(M.base_ring(), M.nrows())
    reduction_points = [pt for pt,p in poincare_map.items() if p >= 1]
    reduction_points.sort()
    if reduction_points == []:
        logger.info("already fuchsian")
    else:
        for pt in reduction_points:
            logger.info("rank = %d, x = %s" % (poincare_map[pt], pt))
    while reduction_points:
        pointidx = rng.randint(0, len(reduction_points) - 1)
        point = reduction_points[pointidx]
        prank = poincare_map[point]
        if prank < 1:
            del reduction_points[pointidx]
            continue
        while True:
            A0 = matrix_c0(M, x, point, prank)
            if A0.is_zero(): break
            A1 = matrix_c1(M, x, point, prank)
            try:
                U, V = alg1x(A0, A1, x)
            except FuchsiaError as e:
                logger.debug("Managed to fuchsify matrix to this state:\n"
                        "%s\nfurther reduction is pointless:\n%s" % (M, e))
                raise FuchsiaError("matrix cannot be reduced to Fuchsian form")
            try:
                point2, P, M = min(iter_reductions(point, U), \
                        key=lambda point2, P, M2: matrix_complexity(M2))
            except ValueError as e:
                point2 = any_integer(rng, M.base_ring(), poincare_map)
                P = fuchsia_simplify(U*V, x)
                M = balance_transform(M, P, point, point2, x)
                M = fuchsia_simplify(M, x)
                logger.info("Will introduce an apparent singularity at %s." % point2)
            logger.debug(
                "Applying balance between %s and %s with projector:\n%s" % (point, point2, P))
            combinedT = combinedT * balance(P, point, point2, x)
            if point2 not in poincare_map:
                poincare_map[point2] = 1
        poincare_map[point] = prank - 1
    combinedT = fuchsia_simplify(combinedT, x)
    logger.exit("fuchsify")
    return M, combinedT

def fuchsify_off_diagonal_blocks(m, x, eps, b=None):
    logger.enter("fuchsify_off_diagonal_blocks")
    n = m.nrows()
    if b is None:
        m, t, b = block_triangular_form(m)
    else:
        t = identity_matrix(SR, n)
    for i, (ki, ni) in enumerate(b):
        for j, (kj, nj) in enumerate(reversed(b[:i])):
            pts = singularities(m.submatrix(ki, kj, ni, nj), x)
            printed = False
            while any(pts.values()):
                bj = m.submatrix(ki, kj, ni, nj)
                if bj.is_zero():
                    break
                for x0, p in pts.items():
                    if p < 1:
                        continue
                    if not printed:
                        printed = True
                        logger.info("Fuchsifying block (%d, %d) (%d, %d)" % (ki, ni, kj, nj))
                        logger.debug("  singular points = %s" % (pts,))
                    a0 = matrix_residue(m.submatrix(ki, ki, ni, ni)/eps, x, x0)
                    b0 = matrix_c0(bj, x, x0, p)
                    c0 = matrix_residue(m.submatrix(kj, kj, nj, nj)/eps, x, x0)

                    d_vars = [gensym() for i in range(ni*nj)]
                    d = matrix(SR, ni, nj, d_vars)
                    eq = d + eps/p*(a0*d - d*c0) + b0/p
                    sol = solve(eq.list(), *d_vars, solution_dict=True)
                    d = d.subs(sol[0])

                    t0 = identity_matrix(SR, n)
                    t0[ki:ki+ni, kj:kj+nj] = \
                            d/(x-x0)**p if not (x0 == oo) else d*(x**p)
                    m = fuchsia_simplify(transform(m, x, t0), x)

                    t = fuchsia_simplify(t*t0, x)
                    pts[x0] -= 1
    logger.exit("fuchsify_off_diagonal_blocks")
    return m, t

def reduce_at_one_point(M, x, v, p, v2=oo):
    """Given a system of differential equations of the form dF/dx=M*F,
    with M having a singularity around x=v with Poincare rank p,
    try to find a transformation T, which will reduce p by 1 (but
    possibly introduce another singularity at x=v2). Return the
    transformed M and T.
    """
    assert M.is_square()
    assert p > 0
    n = M.nrows()
    combinedT = identity_matrix(M.base_ring(), n)
    while True:
        A0 = matrix_c0(M, x, v, p)
        if A0.is_zero(): break
        A1 = matrix_c1(M, x, v, p)
        U, V = alg1x(A0, A1, x)
        P = U*V
        M = balance_transform(M, P, v, v2, x)
        M = fuchsia_simplify(M, x)
        combinedT = combinedT * balance(P, v, v2, x)
    combinedT = fuchsia_simplify(combinedT, x)
    return M, combinedT

def find_dual_basis_spanning_left_invariant_subspace(A, U, rng):
    """Find matrix V, such that it's rows span a left invariant
    subspace of A and form a dual basis with columns of U (that
    is, V*U=I).
    """
    evlist = []
    for eigenval, eigenvects, evmult in A.eigenvectors_left():
        evlist.extend(eigenvects)
    rng.shuffle(evlist)
    W = matrix(evlist)
    try:
        M = solve_left_fixed(W*U, identity_matrix(U.ncols()))
        return M*W
    except ValueError:
        return None

def alg1x(A0, A1, x):
    #if not matrix_is_nilpotent(A0):
    #    raise FuchsiaError("matrix is irreducible (non-nilpotent residue)")
    # We will rely on undocumented behavior of jordan_form() call:
    # we will take it upon faith that it sorts Jordan cells by
    # their size, and puts the largest ones closer to (0, 0).
    A0J, U = A0.jordan_form(transformation=True)
    invU = U.inverse()
    A0J_cs = jordan_cell_sizes(A0J)
    assert all(A0J_cs[i] >= A0J_cs[i+1] for i in range(len(A0J_cs) - 1))
    ncells = len(A0J_cs)
    A0J_css = [sum(A0J_cs[:i]) for i in range(ncells + 1)]
    nsimplecells = sum(1 if s == 1 else 0 for s in A0J_cs)
    u0 = [U[:,A0J_css[i]] for i in range(ncells)]
    v0t = [invU[A0J_css[i+1]-1,:] for i in range(ncells)]
    L0 = matrix([
        [(v0t[k]*(A1)*u0[l])[0,0] for l in range(ncells)]
        for k in range(ncells)
    ])
    L0 = fuchsia_simplify(L0, x)
    L1 = matrix([
        [(v0t[k]*u0[l])[0,0] for l in range(ncells)]
        for k in range(ncells)
    ])
    assert (L1 - diagonal_matrix(
            [0]*(ncells-nsimplecells) + [1]*nsimplecells)).is_zero()
    #zero_rows = [i for i in range(A0.nrows()) if A0J[i,:].is_zero()]
    #zero_cols = [j for j in range(A0.nrows()) if A0J[:,j].is_zero()]
    #assert len(zero_rows) == len(zero_cols) == ncells
    #_L0 = (invU*A1*U)[zero_rows, zero_cols]
    #assert (L0 - _L0).is_zero()
    lam = SR.symbol()
    if not (L0 - lam*L1).determinant().is_zero():
        raise FuchsiaError("matrix is Moser-irreducible")
    S, D = alg1(L0, A0J_cs)
    I_E = identity_matrix(D.base_ring(), A0.nrows())
    for i in range(ncells):
        for j in range(ncells):
            if not D[i,j].is_zero():
                ni = A0J_css[i]
                nj = A0J_css[j]
                for k in range(min(A0J_cs[i], A0J_cs[j])):
                    I_E[ni+k,nj+k] += D[i,j]
    U_t = U*I_E
    invU_t = U_t.inverse()
    return \
        U_t[:, [A0J_css[i] for i in S]], \
        invU_t[[A0J_css[i] for i in S], :]

def alg1(L0, jordan_cellsizes):
    assert L0.is_square()
    ring = L0.base_ring()
    N = L0.nrows()
    I_N = identity_matrix(N)
    S = set()
    D = zero_matrix(ring, N)
    K = sum(1 if s > 1 else 0 for s in jordan_cellsizes)
    while True:
        Lx = copy(L0)
        for i in S:
            Lx[:, i] = 0
            Lx[i, :] = 0
        fi = None
        c = None
        for i in range(0, N):
            if i not in S:
                try:
                    c = solve_right_fixed(Lx[:,0:i], Lx[:,i])
                except ValueError:
                    # No solution found; vectors are independent.
                    continue
                fi = i
                break
        assert fi is not None
        assert c is not None
        D0 = matrix(ring, N)
        invD0 = matrix(ring, N)
        for j in range(fi):
            D0[j, fi] = -c[j, 0]
            invD0[j, fi] = -c[j, 0] \
                    if jordan_cellsizes[j] == jordan_cellsizes[fi] else 0
        L0 = (I_N - invD0)*L0*(I_N + D0)
        D = D + D0 + D*D0
        if fi < K:
            break
        S.add(fi)
    # Check Alg.1 promise
    for j in range(N):
        for k in range(N):
            if (j not in S) and ((k in S) or k == fi):
                assert L0[j,k] == 0
    S.add(fi)
    return S, D

#==================================================================================================
# Step II: Normalize
#==================================================================================================

def is_normalized(M, x, eps):
    """Return True if (a Fuchsian) matrix M is normalized, that
    is all the eigenvalues of it's residues in x lie in [-1/2, 1/2)
    range (in limit eps->0). Return False otherwise.

    Examples:
    >>> x, e = var("x epsilon")
    >>> is_normalized(matrix([[(1+e)/3/x, 0], [0, e/x]]), x, e)
    True
    """
    points = singularities(M, x)
    for x0, p in points.items():
        M0 = matrix_residue(M, x, x0)
        for ev in M0.eigenvalues():
            ev = limit_fixed(ev, eps, 0)
            if not (Rational((-1, 2)) <= ev and ev < Rational((1, 2))):
                return False
    return True

def are_diagonal_blocks_reduced(m, b, x, eps):
    """Return True if diagonal blocks of `m` are normalized, that is all the eigenvalues of theirs
    residues in `x` lie in the range [-1/2, 1/2) in eps->0 limit; return False otherwise. Diagonal
    blocks are defined by the list `b` which corresponds to the equivalent value returned by the
    `block_triangular_form` routine.

    Examples:
    >>> x, e = var("x epsilon")
    >>> are_diagonal_blocks_reduced(matrix([[(1+e)/3/x, 0], [0, e/x]]), [(0,1),(1,1)], x, e)
    True
    """
    for ki, ni in b:
        mi = fuchsia_simplify(m.submatrix(ki, ki, ni, ni), x)
        if not is_normalized(mi, x, eps):
            return False
    return True

def reduce_diagonal_blocks(m, x, eps, b=None, seed=0):
    """Given a lower block-triangular system of differential equations of the form dF/dx=m*F,
    find a transformation that will shift all eigenvalues of all residues of all its diagonal
    blocks into the range [-1/2, 1/2) in eps->0 limit. Return the transformed matrix m and the
    transformation; raise FuchsiaError if such transformation is not found. Diagonal blocks
     are defined by the list `b` which corresponds to the equivalent value returned by the
    `block_triangular_form` routine.
    """
    logger.enter("reduce_diagonal_blocks")
    n = m.nrows()
    if b is None:
        m, t, b = block_triangular_form(m)
    else:
        t = identity_matrix(SR, n)
    for i, (ki, ni) in enumerate(b):
        mi = fuchsia_simplify(m.submatrix(ki, ki, ni, ni), x)
        ti = identity_matrix(SR, ni)
        logger.info("reducing block #%d (%d,%d)" % (i,ki,ni))
        if logger.is_verbose():
            logger.debug("\n%s" % matrix_str(mi, 2))

        mi_fuchs, ti_fuchs = fuchsify(mi, x, seed=seed)
        ti = ti*ti_fuchs

        mi_norm, ti_norm = normalize(mi_fuchs, x, eps, seed=seed)
        ti = ti*ti_norm

        mi_eps, ti_eps = factorize(mi_norm, x, eps, seed=seed)
        ti = ti*ti_eps

        t[ki:ki+ni, ki:ki+ni] = ti
    mt = fuchsia_simplify(transform(m, x, t), x)
    logger.exit("reduce_diagonal_blocks")
    return mt, t


def cmap_add_div(Cmap, C, pi, ki, x0):
    # M += C*(x-pi)^ki/(x-x0)
    if pi == x0:
        Cmap[pi, ki-1] += C
    elif ki >= 0:
        assert pi == 0
        for k in range(0, ki):
            Cmap[pi, k] += C * x0**(ki-1-k)
        Cmap[x0, -1] += C * x0**ki
    else:
        for k in range(ki, 0):
            Cmap[pi, k] -= C * (x0-pi)**(ki-1-k)
        Cmap[x0, -1] += C * (x0-pi)**ki

def cmap_add_mul(Cmap, C, pi, ki, x0):
    # M += C*(x-pi)^ki*(x-x0) = C*(x-pi)^ki*{(x-pi) + (pi-x0)}
    Cmap[pi, ki] += C*(pi-x0)
    if ki == -1:
        Cmap[SR(0), ki+1] += C
    else:
        Cmap[pi, ki+1] += C

def memoize(f):
    @wraps(f)
    def wrapper(self):
        key = f.__name__+"_cache"
        if not hasattr(self, key):
#            print("cache miss: %s" % key)
            setattr(self, key, f(self))
#        else:
#            print("cache hit: %s" % key)
        return getattr(self, key)
    return wrapper

class FuchsianSystem(object):
    @staticmethod
    def from_M(m, x, eps):
        # singular points
        xs = sorted(singularities(m, x).keys())
        # residues in singular points
        rs = {}
        for xi in xs:
            rs[xi] = SystemResidue(matrix_residue(m, x, xi), eps)
        # TODO: check if system's residues have half-integer eigenvalues.
        # If yes, then inform user about that and let us remove resonances only.
        # Though, such a form does not allow to construct an eps-form, it
        # however can be used to find solutions as a power series using
        # semi-analytical methods.
        return FuchsianSystem(rs, x, eps)

    def __init__(self, rs, x, eps):
        self.rs = rs
        self.xs = sorted(rs.keys())
        self.x = x
        self.eps = eps

    def __str__(self):
        #s = "zeroes:\n  %s\n" % (self.xs,)
        s = ""
        a = 0
        for xi, ri in self.rs.items():
            s += "a = %s (x = %s)\n%s\n" % (a, xi, ri)
            a += 1
        return s

    def write(self, filename):
        m = SR(0)
        for xi, ri in self.rs.items():
            if xi == oo:
                continue
            m += ri.m/(self.x-xi)
        export_matrix_to_file(filename, m)

    def _balance(self, i1, j1, i2, j2):
        x1, x2 = self.xs[i1], self.xs[i2]
        r1, r2 = self.rs[x1], self.rs[x2]
        v1, v2 = r1.eigenvectors_right()[j1][1], r2.eigenvectors_left()[j2][1]
        print("x1 = %s" % x1)
        print("x2 = %s" % x2)
        print("v1 = %s" % v1)
        print("v2 = %s" % v2)

        sp = simple(dot_product(v1, v2))
        print("sp = %s" % sp)
        P = simple(cross_product(v1, v2) / sp)
        print("P  =\n%s" % P)
        m = balance_transform(self.m, P, x1, x2, self.x)
        return FuchsianSystem.from_M(m, self.x, self.eps)

    #@logcall
    def balance(self, i1, j1, i2, j2):
        x1, x2 = self.xs[i1], self.xs[i2]
        r1, r2 = self.rs[x1], self.rs[x2]
        v1, v2 = r1.eigenvectors_right()[j1][1], r2.eigenvectors_left()[j2][1]

        sp = simple(dot_product(v1, v2))
        print("sp = %s" % sp)
        P = simple(cross_product(v1, v2) / sp)
        print("P  =\n%s" % P)

        coP = 1-P
        if x1 == oo:
            Cmap = defaultdict(SR.zero)
            for pi, Ci in self.rs.items():
                if pi == oo:
                    continue
                # coP Ci coP (x-pi)^-1 + P Ci P (x-pi)^-1
                Cmap[pi, -1] += coP*Ci.m*coP + P*Ci.m*P
                # coP Ci P -(x-x2) (x-pi)^-1
                cmap_add_mul(Cmap, -coP*Ci.m*P, pi, -1, x2)
                # P Ci coP -1/(x-x2) (x-pi)^-1
                cmap_add_div(Cmap, -P*Ci.m*coP, pi, -1, x2)
            # -P/(x-x2)
            Cmap[x2, -1] -= P
        elif x2 == oo:
            Cmap = defaultdict(SR.zero)
            for pi, Ci in self.rs.items():
                if pi == oo:
                    continue
                # coP Ci coP (x-pi)^-1 + P Ci P (x-pi)^-1
                Cmap[pi, -1] += coP*Ci.m*coP + P*Ci.m*P
                # P Ci coP -(x-x1) (x-pi)^-1
                cmap_add_mul(Cmap, -P*Ci.m*coP, pi, -1, x1)
                # coP Ci P -1/(x-x1) (x-pi)^-1
                cmap_add_div(Cmap, -coP*Ci.m*P, pi, -1, x1)
            # P/(x-x1)
            Cmap[x1, -1] += P
        else:
            Cmap = defaultdict(SR.zero, [((pi,-1),ri.m) for pi,ri in self.rs.items()])
            for pi, Ci in self.rs.items():
                if pi == oo:
                    continue
                # coP Ci P (x1-x2)/(x-x1) (x-pi)^-1
                cmap_add_div(Cmap, coP*Ci.m*P*(x1-x2), pi, -1, x1)
                # P Ci coP (x2-x1)/(x-x2) (x-pi)^-1
                cmap_add_div(Cmap, P*Ci.m*coP*(x2-x1), pi, -1, x2)
            # P/(x-x1) - P/(x-x2)
            Cmap[x1, -1] += P
            Cmap[x2, -1] -= P
        rs = {}
        Coo = SR.zero()
        for key, C in Cmap.items():
            pi, ki = key
            if pi == oo:
                continue
            Cmap[key] = C = simple(C)
            if C.is_zero():
                del Cmap[key]
            else:
                assert(ki == -1)
                rs[pi] = SystemResidue(C)
                Coo -= C
        Coo = simple(Coo)
        rs[oo] = SystemResidue(Coo)
        #T = fuchsia_simplify(self.T * balance(P, x1, x2, self.x))
        #return RationalSystem(dict(Cmap), self.x, T, self.trace + [("balance", P, x1, x2)])
        return FuchsianSystem(rs, self.x, self.eps) 

    #@logcall
    @memoize
    def balances(self):
        res = []
        for i1,x1 in enumerate(self.xs):
            for i2,x2 in enumerate(self.xs):
                if i1 == i2:
                    continue
                for j1, j2, aux in self.balances2(x1, x2):
                    res.append((i1,j1, i2,j2, aux))
        return res

    #@logcall
    def balances2(self, x1, x2):
        #assert(not (x1-x2).is_zero())
        res = []
        r1, r2 = self.rs[x1], self.rs[x2]
        for j1, (a1,v1) in enumerate(r1.eigenvectors_right()):
            for j2, (a2,v2) in enumerate(r2.eigenvectors_left()):
                sp = simple(dot_product(v1, v2))
                if sp.is_zero():
                    continue
                res.append((j1, j2, (a1,a2, sp)))
        return res

    #@logcall
    @memoize
    def eigenvalues(self):
        res = []
        for xi in self.xs:
            res.append((xi, self.rs[xi].eigenvalues()))
        return res

#    @memoize
#    @logcall
#    def eigenvectors_left(self):
#        return self._eigenvectors(left=True)
#
#    @memoize
#    @logcall
#    def eigenvectors_right(self):
#        return self._eigenvectors(left=False)
#
#    def _eigenvectors(self, left=True):
#        res = []
#        for xi, ri in self.rs.items():
#            eigenvectors = ri.eigenvectors_left
#            if left is not True:
#                eigenvectors = ri.eigenvectors_right
#            res.append((xi, eigenvectors()))
#        return res

    def normalize(self):
        pass

class SystemResidue(object):
    def __init__(self, m, eps=None):
        self.m = m
        self.eps = eps

    def __str__(self):
        return str(self.m)

    def __repr__(self):
        return self.m.__repr__()

    @logcall
    @memoize
    def eigenvalues(self):
        return simple(mathematica.Eigenvalues(self.m).sage())

    @logcall
    @memoize
    def eigenvectors_left(self):
        return self._eigenvectors(left=True)

    @logcall
    @memoize
    def eigenvectors_right(self):
        return self._eigenvectors(left=False)

    def _eigenvectors(self, left=True):
        m = self.m.transpose() if left is True else self.m
        s = simple(mathematica.Eigensystem(m).sage())
        res = zip(s[0], s[1])
        return res

#    def _eigenvectors(self, left=True):
#        res = []
#        eigenvectors = self.m.eigenvectors_left
#        if left is not True:
#            eigenvectors = self.m.eigenvectors_right
#        res = [vvm[0:2] for vvm in eigenvectors()]
#        return res

# helpers
@logcall
def simple(obj, x=None):
    def _simplify(ex):
            return mathematica.Simplify(ex).sage()
    if hasattr(obj, "apply_map"):
        return obj.apply_map(_simplify)
    else:
        return _simplify(obj)

def ex_in(ex, lst):
    #if any((ex-t).is_zero() for t in printed):
    for t in lst:
        if (ex-t).is_zero():
            return True
    return False
        

class NormalizeAssistant(object):
    def __init__(self, m):
        self.m = m           # current matrix to balance
        self.history = []    # matrices balanced so far

    # Entry point
    def start(self):
#        cmd = self.cmd_main
#        while cmd != None:
#            cmd = cmd()
#        return
        try:
            cmd = self.cmd_main
            while cmd != None:
                cmd = cmd()
        except Exception as error:
            print("UNEXPECTED ERROR:\n%s" % (error,))
            return self.start()

    # Helpers
    def choose_balance(self, i1, i2):
        x1, x2 = self.m.xs[i1], self.m.xs[i2]
        bs = self.m.balances2(x1, x2)
        for i, (j1, j2, aux) in enumerate(bs):
            print("[%s] %s" % (i, aux))
        n = raw_input("(choose balance) > ")
        try:
            n = int(n)
        except ValueError:
            return None
        if n < 0 or n > len(bs):
            return None
        b = bs[n]
        return (i1, b[0], i2, b[1])

    # Commands
    def cmd_balance(self, b):
        print("balancing with %s" % (b,))

        mm = apply(self.m.balance, b)
        self.history.append(mm)
        self.m = mm
        return self.cmd_main

    def cmd_balance_res(self, args):
        if len(args) != 2:
            return self.cmd_error("needs exactly 2 arguments")
        try:
            i1, i2 = [int(n) for n in args]
        except ValueError:
            return self.cmd_error("arguments out of range (type e for valid values)")
        b = self.choose_balance(i1, i2)
        if b is None:
            return self.cmd_error("unknown balance")
        return self.cmd_balance(b[0:4])

    def cmd_error(self, msg):
        print("ERROR: %s" % msg)
        return self.cmd_main

    def cmd_main(self):
        cmd = raw_input("> ").split(' ')
        cmd, args = cmd[0], cmd[1:]
        if cmd == 's':
            return self.cmd_save(args)
        elif cmd == 'b':
            return self.cmd_balance_res(args)
        elif cmd == 'e':
            return self.cmd_print_eigenvalues
        elif cmd == 'h':
            return self.cmd_print_help
        elif cmd == 'q':
            return None
        elif cmd == '':
            return self.cmd_print_help
        return self.cmd_unknown_command(cmd)

    def cmd_print_eigenvalues(self):
        for i, (xi, a) in enumerate(self.m.eigenvalues()):
            print("[%s] x = %s\n    %s" % (i, xi, a))
        return self.cmd_main

    def cmd_print_help(self):
        """
        Available commands:
            b <i> <j>    balance xi, xj points
            e            show eigenvectors
            h            print this help
            q            quit
            s <file>     save matrix to file
        """
        print(inspect.cleandoc(self.cmd_print_help.__doc__))
        return self.cmd_main

    def cmd_save(self, args):
        if len(args) != 1:
            return self.cmd_error("needs exactly 1 arguments")
        name = args[0]
        self.m.write(name)
        return self.cmd_main

    def cmd_unknown_balance(self, n):
        print("Unknown balance: '%s'" % n)
        return self.cmd_main

    def cmd_unknown_command(self, cmd):
        print("Unknown command: '%s'" % cmd)
        return self.cmd_print_help

def normalize(m, x, eps, seed=0):
    """Given a Fuchsian system of differential equations of the
    form dF/dx=m*F, find a transformation that will shift all
    the eigenvalues of m's residues into [-1/2, 1/2) range (in
    the limit eps->0). Return the transformed matrix m and the
    transformation. Raise FuchsiaError if such transformation
    is not found.
    """
    class State(object):
        def __init__(self, m, x, eps, seed):
            # random
            self.random = Random(seed)
            # points
            self.points = singularities(m, x).keys()
            # ev_cum
            self.ev_cum = {}
            for x0 in self.points:
                pos, neg = 0, 0
                a0 = matrix_residue(m, x, x0)
                for ev in a0.eigenvalues():
                    ev0 = limit_fixed(ev, eps, 0)
                    if ev0 > 0:
                        pos += ev0
                    if ev0 < 0:
                        neg += ev0
                    logger.debug("x = %s, ev = %s" % (x0, ev0))
                self.ev_cum[x0] = [pos,neg]
            # x0
            logger.debug("state.ev_cum = %s" % self.ev_cum)
            self.x0 = None

        def is_normalized(self):
            for ev in self.ev_cum.values():
                if ev != [0,0]:
                    return False
            return True

        def pairs(self):
            points = [x0 for x0 in self.points if self.ev_cum[x0] != [0,0]]
            if len(points) == 1:
                self.select_x0(points[0])
                points.append(self.x0)
            return permutations(points, 2)

        def select_x0(self, x1):
            if self.x0 is not None:
                return
            for x0 in self.points:
                if x0 != x1:
                    self.x0 = x0
                    return

        def update_ev_cum(self, x2, x1):
            ev_cum_x1, ev_cum_x2 = self.ev_cum[x1], self.ev_cum[x2]
            if ev_cum_x1[0] > 0:
                ev_cum_x1[0] -= 1
            else:
                ev_cum_x1[1] -= 1
            if ev_cum_x2[1] < 0:
                ev_cum_x2[1] += 1
            else:
                ev_cum_x2[0] += 1

    logger.enter("normalize")
    state = State(m, x, eps, seed)
    T = identity_matrix(m.base_ring(), m.nrows())
    if state.is_normalized():
        logger.info("already normalized")
    i = 0
    while not state.is_normalized():
        i += 1
        m = fuchsia_simplify(m, x)
        logger.info("step %s" % i)
        balances = find_balances(m, x, eps, state)
        b = select_balance(balances, eps, state)
        if b is None:
            raise FuchsiaError("can not balance matrix")
        logger.info("  balancing x = %s and x = %s" % (b[1],b[2]))
        if logger.is_verbose():
            logger.debug("\n      use the balance:\n        %s\n" % b)

        cond, x1, x2, a0_eval, b0_eval, a0_evec, b0_evec, scale = b
        if cond == 1:
            P = cross_product(a0_evec, b0_evec) / scale
            m = balance_transform(m, P, x1, x2, x)
            T0 = balance(P, x1, x2, x)
            state.update_ev_cum(x1, x2)
        else:
            P = cross_product(b0_evec, a0_evec) / scale
            m = balance_transform(m, P, x2, x1, x)
            T0 = balance(P, x2, x1, x)
            state.update_ev_cum(x2, x1)

        T = fuchsia_simplify(T*T0, x)
    logger.exit("normalize")
    return m, T

def find_balances(m, x, eps, state={}):
    residues = {}
    for x1, x2 in state.pairs():
        logger.debug("trying to balance x = %s and x = %s" % (x1,x2))
        for xi in [x1,x2]:
            if xi not in residues:
                residues[xi] = matrix_residue(m, x, xi)
        a0, b0 = residues[x1], residues[x2]

        a0_evr, b0_evl = eigenvectors_right(a0), eigenvectors_left(b0)

        if logger.is_verbose():
            msg = "\n  Eigenvalues:\n"
            msg += "    x = %s:\n" % x1
            a0_evals = [];
            for ev, evec, emult in a0_evr:
                a0_evals += [ev]*emult
            msg += "        %s\n" % str(a0_evals)
            msg += "    x = %s:\n" % x2
            b0_evals = [];
            for ev, evec, emult in b0_evl:
                b0_evals += [ev]*emult
            msg += "        %s\n" % str(b0_evals)
            logger.debug(msg)

        balances_1 = find_balances_by_cond(a0_evr, b0_evl, lambda a0_eval, b0_eval: limit_fixed(a0_eval, eps, 0) < -0.5)
        for balance in balances_1:
            balance = [1, x1, x2] + balance
            yield balance

        a0_evl, b0_evr = eigenvectors_left(a0), eigenvectors_right(b0)
        balances_2 = find_balances_by_cond(a0_evl, b0_evr, lambda a0_eval, b0_eval: limit_fixed(a0_eval, eps, 0) >= 0.5)
        for balance in balances_2:
            balance = [2, x1, x2] + balance
            yield balance

def find_balances_by_cond(a0_ev, b0_ev, cond):
    res = []
    for a0_eval, a0_evecs, a0_evmult in a0_ev:
        for b0_eval, b0_evecs, b0_evmult in b0_ev:
            if not cond(a0_eval, b0_eval):
                logger.debug("Balance rejected:\n    a0_eval = %s\n    b0_eval = %s" % (a0_eval, b0_eval))
                continue
            for a0_evec in a0_evecs:
                for b0_evec in b0_evecs:
                    scale = fuchsia_simplify(dot_product(a0_evec, b0_evec))
                    balance = [a0_eval, b0_eval, a0_evec, b0_evec, scale]
                    if scale == 0:
                        logger.debug("Balance rejected:\n    a0_eval = %s\n    b0_eval = %s\n    a0_evec = %s\n    b0_evec = %s\n    scale   = %s" % tuple(balance))
                        continue
                    logger.debug("Balance found:\n    a0_eval = %s\n    b0_eval = %s\n    a0_evec = %s\n    b0_evec = %s\n    scale   = %s" % tuple(balance))
                    res.append(balance)
    return res

def select_balance(balances, eps, state={}):
    min_degree, min_balance = None, None
    bs = []
    for b in balances:
        cond, x1, x2, a0_eval, b0_eval, a0_evec, b0_evec, scale = b
        if (cond == 1) and limit_fixed(a0_eval, eps, 0) < -0.5 and \
                limit_fixed(b0_eval, eps, 0) >= 0.5:
            degree = max(scale.numerator().degree(eps), scale.denominator().degree(eps))
            if degree < 4:
                return b
            if (min_degree is None) or (min_degree > degree):
                 min_degree = degree
                 min_balance = b
        elif (cond == 2) and limit_fixed(a0_eval, eps, 0) >= 0.5 and \
                limit_fixed(b0_eval, eps, 0) < -0.5:
            degree = max(scale.numerator().degree(eps), scale.denominator().degree(eps))
            if degree < 4:
                return b
            if (min_degree is None) or (min_degree > degree):
                 min_degree = degree
                 min_balance = b
        bs.append(b)
    if min_balance is not None:
        return min_balance

    x0 = state.x0
    if x0 is None:
        for b in bs:
            cond, x1, x2, ev1, ev2 = b[:5]
            if cond == 1:
                x0 = x2
                break
            if cond == 2:
                x0 = x1
                break
        logger.info("  select x0 = %s" % x0)
        state.x0 = x0

    balances_x0 = [b for b in bs if (b[0] == 1 and b[2] == x0) or (b[0] == 2 and b[1] == x0)]
    b = state.random.choice(balances_x0) if balances_x0 else None
    return b

def eigenvectors_left(m):
    if m._cache is None:
        m._cache = {}
    if "eigenvectors_left" not in m._cache:
        res = simplify(m.eigenvectors_left())
        m._cache["eigenvectors_left"] = res
    return m._cache["eigenvectors_left"]

def eigenvectors_right(m):
    if m._cache is None:
        m._cache = {}
    if "eigenvectors_right" not in m._cache:
        res = m.eigenvectors_right()
        m._cache["eigenvectors_right"] = res
    return m._cache["eigenvectors_right"]

#==================================================================================================
# Step III: Factorize
#==================================================================================================

def gensym():
    sym = SR.symbol()
    SR.symbols[str(sym)] = sym
    return sym

def factorize(M, x, epsilon, b=None, seed=0):
    """Given a normalized Fuchsian system of differential equations:
        dF/dx = M(x,epsilon)*F,
    try to find a transformation that will factor out an epsilon
    from M. Return a transformed M (proportional to epsilon)
    and T. Raise FuchsiaError if epsilon can not be factored.
    """
    logger.info("-> factorize")
    n = M.nrows()
    M = fuchsia_simplify(M, x)
    if epsilon not in expand(M/epsilon).variables():
        logger.info("   already in epsilon form")
        logger.info("<- factorize")
        return M, identity_matrix(SR, n)
    rng = Random(seed)
    mu = gensym()
    if b is None:
        T_symbols = [gensym() for i in range(n*n)]
        T = matrix(SR, n, n, T_symbols)
    else:
        T, T_symbols = matrix(SR, n), []
        for ki,ni in b:
            for i in range(ki,ki+ni):
                for j in range(ki+ni):
                    sym = gensym()
                    T[i,j] = sym
                    T_symbols.append(sym)
    eqs = []
    for point, prank in singularities(M, x).items():
        assert prank == 0
        logger.debug("    processing point x = %s" % point)
        R = matrix_c0(M, x, point, 0)
        R = fuchsia_simplify(R)
        eq = (R/epsilon)*T-T*(R.subs({epsilon: mu})/mu)
        eq = fuchsia_simplify(eq)
        eqs.extend(eq.list())
    logger.info("    found %d equations with %d unknowns" % (len(eqs), len(T_symbols)))
    solutions = fuchsia_solve(eqs, T_symbols)
    for solution in solutions:
        S = T.subs(solution)
        # Right now S likely has a number of free variables in
        # it; we can set them to arbibtrary values, as long as
        # it'll make S invertible.
        rndrange = 0
        while True:
            try:
                sT = S.subs([
                    e==rng.randint(-rndrange, rndrange)
                    for e in S.variables() if e != epsilon
                ])
                sT = fuchsia_simplify(sT,x)
                M = fuchsia_simplify(transform(M, x, sT), x)
                # We're leaking a bunch of temprary variables here,
                # which accumulate in SR.variables, but who cares?
                logger.info("<- factorize")
                return M, sT
            except (ZeroDivisionError, ValueError):
                rndrange += 1 + rndrange//16
                # We've tried a bunch of substituions, and they didn't
                # work. Is the matrix at all invertible? Let's check.
                if rndrange == 16 and not S.is_invertible():
                    break
    raise FuchsiaError("can not factor epsilon")

#==================================================================================================
# Helpers
#==================================================================================================

def matrix_complexity(M):
    return len(str(M.list()))

def simplify_by_jordanification(M, x):
    """Try to simplify matrix M by constant transformations that
    transform M's residues into their Jordan forms. Return the
    simplified matrix and the transformation. If none of the
    attempted transformations reduce M's complexity (as measured
    by 'matrix_complexity()'), return the original matrix and
    the identity transformation.
    """
    minM = M
    minC = matrix_complexity(M)
    minT = identity_matrix(M.base_ring(), M.nrows())
    for point, prank in singularities(M, x).items():
        R = matrix_c0(M, x, point, prank)
        J, T = R.jordan_form(transformation=True)
        MM = fuchsia_simplify(transform(M, x, T), x)
        C = matrix_complexity(MM)
        if C < minC:
            minM = MM
            minC = C
            minT = T
    return minM, minT

def common_factor(expressions, filter):
    """Factorize given expressions, select those factors for
    which 'filter(factor)' is True, and return the product of
    factors common to all the expressions.

    Examples:
    >>> x = var("x")
    >>> common_factor([x*x-1, x+1], lambda f: True)
    x + 1
    >>> common_factor([1/x**2, 2/x**3, 3/x**4], lambda f: True)
    x^(-2)

    Note that if there is a mix of positive and negative exponents
    of a given factor, this function will use (one of) the most
    frequently occurring exponent:
    >>> common_factor([x, 1/x, 2/x**2, 3/x], lambda f: True)
    1/x
    """
    factor2exp2count = defaultdict(lambda: defaultdict(lambda: 0))
    for i, expr in enumerate(expressions):
        factors = dict(expr.factor_list())
        for factor, n in factors.items():
            if not filter(factor): continue
            if factor in factor2exp2count:
                factor2exp2count[factor][n] += 1
            else:
                if i > 0: factor2exp2count[factor][0] = i
                factor2exp2count[factor][n] = 1
        for factor, exps in factor2exp2count.items():
            if factor not in factors:
                exps[0] += 1
    result = SR(1)
    for factor, exp2count in factor2exp2count.items():
        exps = exp2count.keys()
        minn = min(exps)
        maxn = max(exps)
        if minn > 0: result *= factor**minn
        if maxn < 0: result *= factor**maxn
        if minn <= 0 and maxn >= 0:
            bestn = max(exps, key=lambda exp: exp2count[exp])
            result *= factor**bestn
    return result

def simplify_by_factorization(M, x):
    """Try to simplify matrix M by a constant transformation
    that extracts common factors found in M (if any). Return
    the simplified matrix and the transformation.
    """
    n = M.nrows()
    T = identity_matrix(M.base_ring(), n)
    factors = []
    for i in range(n):
        factor = common_factor(
            [M[i,k] for k in range(n) if i != k and not M[i,k].is_zero()] +
            [1/M[k,i] for k in range(n) if k != i and not M[k,i].is_zero()],
            lambda e: x not in e.variables())
        if factor != 1:
            dT = identity_matrix(M.base_ring(), n)
            dT[i,i] = T[i,i] = factor
            M = transform(M, x, dT)
    logger.debug(
            "stripping common factors with this transform:\ndiagonal_matrix([\n  %s\n])" % ",\n  ".join(str(e) for e in T.diagonal()))
    return fuchsia_simplify(M, x), T

#==============================================================================
# Import/Export routines
#==============================================================================

_parser = Parser(make_int=SR, make_float=SR,
        make_var=lambda v: I if v in "Ii" else SR.var(v),
        make_function={"sqrt": sqrt})

def parse(s):
    """Convert a given string to a Sage expression.

    >>> parse("(1 + I*eps)*(1 - I*eps)").simplify_rational()
    eps^2 + 1
    """
    return _parser.parse(s)

def import_matrix_from_file(filename):
    """Read and return a matrix from a named file. Both Mathematica
    and MatrixMarket array formats are supported. The exact format
    will be autodetected.
    """
    with open(filename, 'r') as f:
        hd = f.read(2)
        f.seek(0)
        if hd == "%%":
            return import_matrix_matrixmarket(f)
        elif hd == "{{":
            return import_matrix_mathematica(f)
        else:
            raise ValueError("File '%s' is not in a known matrix format" % filename)

def import_matrix_mathematica(f):
    """Read and return a matrix, stored in the Mathematica format,
    from a file-like object.
    """
    data =  f.read()
    data = data.translate(None, ' \n')
    nr = data.count('},{')+1
    data = data.translate(None, '{}')
    data = data.replace(',', '\n')
    data = data.replace('Sqrt[', 'sqrt(').replace(']', ')')
    nc = int((data.count('\n')+1)/nr)

    import StringIO
    sio = StringIO.StringIO("%d %d\n" % (nc, nr) + data)
    return import_matrix_matrixmarket(sio).T

def import_matrix_matrixmarket(f):
    """Read and return a matrix, stored in the MatrixMarket
    array format, from a file-like object.
    """
    lineno = 0
    while True:
        s = f.readline()
        lineno += 1
        if not s.startswith('%'):
            break
    try:
        nrows, ncols = map(int, s.split())
    except ValueError:
        msg = ("Can not determine size of the matrix\n"
               "Make sure the file is in the MatrixMarket array format")
        raise FuchsiaError(msg)
    data = []
    for s in f.readlines():
        lineno += 1
        if s.startswith('%'):
            continue
        try:
            expr = parse(s)
        except SyntaxError:
            msg = ("Malformad expression at line %d:\n%s\n"
                   "Make sure the file is in the MatrixMarket array format")
            raise FuchsiaError(msg % (lineno, s))
        data.append(expr)
    m = matrix(ncols, nrows, data).T
    return m

def export_matrix_to_file(filename, m, fmt="m"):
    """Save matrix into a file in a given matrix format. Set
    'fmt' to "mtx" for MatrixMarket array format, or to "m"
    for Mathematica format.
    """
    with open(filename, 'w') as f:
        if fmt == "mtx":
            export_matrix_matrixmarket(f, m)
        elif fmt == "m":
            export_matrix_mathematica(f, m)
        else:
            raise FuchsiaError("Unknown matrix format '%s'" % fmt)

def export_matrix_mathematica(f, m):
    """Write matrix to a file-like object using Mathematica format."""
    i = 0
    f.write("{")
    nc, nr = m.ncols(), m.nrows()
    for row in m.rows():
        i += 1
        f.write("{")
        j = 0
        for mij in row:
            j += 1
            f.write(str(mij).replace(' ', ''))
            if j < nc:
                f.write(',')
        f.write("}")
        if i < nr:
            f.write(",")
    f.write('}')

def export_matrix_matrixmarket(f, m):
    """Write matrix to a file-like object using MatrixMarket
    array format.
    """
    f.write("%%MatrixMarket matrix array Fuchsia[symbolic] general\n")
    f.write("%d %d\n" % (m.nrows(), m.ncols()))
    i = 0
    for col in m.columns():
        i += 1
        f.write('%% #%d\n' % i)
        for mij in col:
            f.write(str(mij).replace(' ', ''))
            f.write('\n')

#==================================================================================================

def main():
    import getopt
    from   contextlib import contextmanager

    @contextmanager
    def profile(file_stats):
        if file_stats is None:
            yield
            return

        import cProfile, pstats
        profile = cProfile.Profile()
        profile.enable()
        try:
            yield
        finally:
            profile.disable()
            with open(file_stats, 'w') as f:
                stats = pstats.Stats(profile, stream=f).sort_stats("cumulative")
                stats.print_stats(50)

    def usage():
        print(__doc__)
        exit(0)

    try:
        print('\033[35;1mFuchsia v%s\033[0m' % (__version__))
        print("Authors: %s\n" % __author__)
        mpath = tpath = profpath = fmt = None
        M = T = None
        x, y, epsilon = SR.var("x y eps")
        fmt = "m"
        seed = 0
        logger.setLevel(logging.INFO)
        kwargs, args = getopt.gnu_getopt(sys.argv[1:],
                "hvl:f:P:x:y:e:m:t:s:", ["help", "use-maple"])
        for key, value in kwargs:
            if key in ["-h", "--help"]: usage()
            if key == "-f":
                if value not in ["m", "mtx"]:
                    raise getopt.GetoptError(
                            "'%s' is not a known matrix format; "
                            "'m' and 'mtx' are supported" % value)
                fmt = value
            if key == "-l":
                fh = logging.FileHandler(value, "w")
                logger_format = '%(levelname)-8s [%(asctime)s] %(message)s'
                fh.setFormatter(ElapsedTimeFormatter(logger_format))
                logger.addHandler(fh)
            if key == "-v": logger.setLevel(logging.DEBUG)
            if key == "-P": profpath = value
            if key == "-x": x = parse(value)
            if key == "-y": y = parse(value)
            if key == "-e": epsilon = SR.var(value)
            if key == "-m": mpath = value
            if key == "-t": tpath = value
            if key == "-s": seed = int(value)
            if key == "--use-maple":
                global USE_MAPLE
                USE_MAPLE = True
        with profile(profpath):
            if len(args) == 2 and args[0] == 'fuchsify':
                M = import_matrix_from_file(args[1])
                M, t1 = simplify_by_factorization(M, x)
                M, t2 = fuchsify(M, x, seed=seed)
                T = t1*t2
            elif len(args) == 2 and args[0] == 'normalize':
                M = import_matrix_from_file(args[1])
                M, t1 = simplify_by_factorization(M, x)
                M, t2 = normalize(M, x, epsilon, seed=seed)
                T = t1*t2
            elif len(args) == 2 and args[0] == 'factorize':
                M = import_matrix_from_file(args[1])
                M, T = factorize(M, x, epsilon, seed=seed)
            elif len(args) == 2 and args[0] == 'sort':
                M = import_matrix_from_file(args[1])
                M, T, B = block_triangular_form(M)
            elif len(args) == 2 and args[0] == 'reduce':
                m = import_matrix_from_file(args[1])
                M, T = epsilon_form(m, x, epsilon, seed=seed)
            elif len(args) == 2 and args[0] == 'cat':
                M = import_matrix_from_file(args[1])
                T = None
            elif len(args) == 3 and args[0] == 'transform':
                M = import_matrix_from_file(args[1])
                t = import_matrix_from_file(args[2])
                M = transform(M, x, t)
            elif len(args) == 3 and args[0] == 'changevar':
                M = import_matrix_from_file(args[1])
                fy = parse(args[2])
                extravars = set(fy.variables()) - set(M.variables())
                if y not in extravars:
                    raise getopt.GetoptError(
                            "'%s' is missing in '%s'" % (y, fy))
                M = change_variable(M, x, y, fy)
                x = y
            elif len(args) == 2 and args[0] == 'info':
                M = import_matrix_from_file(args[1])
                print("Matrix size: %s" % M.nrows())
                print("Matrix complexity: %s" % matrix_complexity(M))
                print("Matrix expansion:")
                for point, prank in singularities(M, x).items():
                    c0 = matrix_c0(M, x, point, prank)
                    print("  C0[%s=%s; %s] = {rank=%s, complexity=%s}" % (
                        x, point, prank, c0.rank(), matrix_complexity(c0)
                    ))
                    if logger.isEnabledFor(logging.DEBUG):
                        for eigenval in sorted(set(c0.eigenvalues()), key=str):
                            print("    eigenvalue: %s" % (eigenval))
                M = None
            elif len(args) == 0:
                usage()
            else:
                raise getopt.GetoptError("unknown command: %s" % (args))
        if M is not None:
            M = partial_fraction(M, x)
            if mpath is not None:
                export_matrix_to_file(mpath, M, fmt=fmt)
            else:
                print(M)
        if tpath is not None and T is not None:
            T = partial_fraction(T, x)
            export_matrix_to_file(tpath, T, fmt=fmt)
    except getopt.GetoptError as error:
        logger.error("%s", error)
        exit(1)
    except FuchsiaError as e:
        msg, tab = str(e), '  '
        msg = tab+msg.replace('\n', '\n'+tab)
        logger.error(msg)
        exit(1)

if __name__ == '__main__':
    main()
