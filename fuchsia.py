#!/usr/bin/env sage
"""\
Usage:
    fuchsia [-hv] [-f <fmt>] [-l <path>] [-P <path>] <command> <args>...

Commands:
    fuchsify [-x <name>] [-m <path>] [-t <path>] <matrix>
        find a transformation that will transform a given matrix
        into Fuchsian form

    normalize [-x <name>] [-e <name>] [-m <path>] [-t <path>] <matrix>
        find a transformation that will transform a given Fuchsian
        matrix into normalized form

    factorize [-x <name>] [-e <name>] [-m <path>] [-t <path>] <matrix>
        find a transformation that will make a given normalized
        matrix proportional to the infinitesimal parameter

    sort [-m <path>] [-t <path>] <matrix>
        find a block-triangular form of the given matrix

    reduce [-x <name>] [-e <name>] [-m <path>] [-t <path>] <matrix>
        find a canonical form of the given matrix

    transform [-x <name>] [-m <path>] <matrix> <transform>
        transform a given matrix using a given transformation

    transform [-m <path>] -y <expr> <matrix>
        change a free variable according to the given expression

Options:
    -h          show this help message
    -f <fmt>    matrix file format: mtx or m (default: mtx)
    -l <path>   save log information to this file
    -v          be verbose and print log information
    -P <path>   save profile report into this file
    -x <name>   use this name for the free variable (default: x)
    -y <expr>   expression for the new free variable (example: "x=(1-y^2)/(1+y^2)")
    -e <name>   use this name for the infinitesimal parameter (default: eps)
    -m <path>   save the resulting matrix into this file
    -t <path>   save the resulting transformation into this file

Arguments:
    <matrix>    read the input matrix from this file
    <transform> read the transformation matrix from this file
"""
from   collections import defaultdict
from   itertools import combinations
import logging
import os.path as path
from   random import Random

from   sage.all import *
from   sage.misc.parser import Parser
from   sage.libs.ecl import ecl_eval
ecl_eval("(ext:set-limit 'ext:heap-size 0)")

logging.basicConfig(
    format='\033[32m%(levelname)s [%(asctime)s]\033[0m\n%(message)s',
    datefmt='%Y-%m-%d %I:%M:%S',
)
logger = logging.getLogger('fuchsia')
logger_format = '%(levelname)s [%(asctime)s]\n%(message)s'

def is_verbose():
    return logger.isEnabledFor(logging.INFO)

__author__ = "Oleksandr Gituliar, Vitaly Magerya"
__author_email__ = "fuchsia@gituliar.org"
__version__ = open(path.join(path.dirname(__file__),"VERSION")).readline().rstrip('\n')
try:
    __commit__ = open(path.join(path.dirname(__file__),"COMMIT")).readline().rstrip('\n')
except IOError:
    __commit__ = "unknown"

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

def limit_fixed(expr, x, x0, cas="maxima"):
    if cas == "maxima":
        return limit_fixed_maxima(expr, x, x0)
    elif cas == "maple":
        return limit_fixed_maple(expr, x, x0)
    else:
        raise NotImplemented

def limit_fixed_maple(expr, x, x0):
    res = maple.limit(simplify(expr), x=x0)
    try:
        res = _parser.parse(str(res))
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

def singularities(m, x, cas="maxima"):
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
    m = m.simplify_rational()
    result = {}
    for expr in m.list():
        points_expr = singularities_expr(expr, x, cas=cas)
        for x0, p in points_expr.iteritems():
            if x0 in result:
                result[x0] = max(result[x0], p)
            else:
                result[x0] = p
    return result

def singularities_expr(expr, x, cas="maxima"):
    if cas == "maple":
        res = singularities_expr_maple(expr, x)
    elif cas == "maxima":
        res = singularities_expr_maxima(expr, x)
    else:
        raise NotImplemented
    return res

def singularities_expr_maple(expr, x):
    if bool(expr == 0):
        return {}

    candidates = {}
    points = [x0 for x0 in f_solve(1/expr, x, cas="maple")]
    for x0 in points:
        if x0 not in candidates:
            candidates[x0] = 0
        else:
            candidates[x0] += 1

    result = {}
    for x0, p in candidates.iteritems():
        while (p >= 0) and bool(limit_fixed(expr*(x-x0)**(p+1), x, x0,  cas="maple") == 0):
            print "REJECTING pole p = %s x = %s" % (p, x0)
            p -= 1
        if p >= 0:
            result[x0] = p

    points = f_solve((1/(expr.subs({x: 1/x})/x**2)).simplify_rational(), x,"maple")
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

def matrix_taylor0(m, x, x0, exp, cas="maxima"):
    """Return the 0-th coefficient of Taylor expansion of
    a matrix M around a finite point x=point, assuming that
    M(x->point)~1/(x-point)**exp.

    Example:
    >>> x = var('x')
    >>> matrix_taylor0(matrix([[x/(x-1), 1/x, x, 1]]), x, 1, 1)
    [1 0 0 0]
    """
    return m.apply_map(lambda ex: limit_fixed(ex*(x-x0)**exp, x, x0, cas))

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

def matrix_c0(M, x, point, p, cas="maxima"):
    """Return the 0-th coefficient of M's expansion at x=point,
    assuming Poincare rank of M at that point is p. If point is
    +Infinity, return the coefficient at the highest power of x.

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
    [0 0]
    [1 0]
    >>> matrix_c0(m*x, x, oo, 2)
    [0 0]
    [1 0]
    """
    if point == oo:
        return matrix_taylor0(M.subs({x: 1/x}), x, 0, p-1, cas)
    else:
        return matrix_taylor0(M, x, point, p+1, cas)

def matrix_c1(M, x, point, p):
    """Return the 1-st coefficient of M's expansion at x=point,
    assuming Poincare rank of M at that point is p. If point is
    +Infinity, return the coefficient at the second-to-highest
    power of x.

    Examples:
    >>> x = var("x")
    >>> m = matrix([[1/x, 1/x**2], [1, 1/(x-1)]])
    >>> matrix_c1(m, x, 0, 1)
    [1 0]
    [0 0]
    >>> matrix_c1(m, x, oo, 1)
    [1 0]
    [0 1]
    """
    if point == oo:
        return matrix_taylor1(M.subs({x: 1/x}), x, 0, p-1)
    else:
        return matrix_taylor1(M, x, point, p+1)

def matrix_residue(M, x, x0, cas="maxima"):
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
    if M._cache.has_key(key):
        return M._cache[key]

    m0 = matrix_c0(M, x, x0, 0, cas)
    if x0 == oo:
        res = -m0
    else:
        res = m0
    M._cache[key] = res
    return res

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
    for i in xrange(J.nrows() - 1):
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
    if logger.isEnabledFor(logging.INFO):
        logger.info("Matrix mask before transformation:\n%s\n" % matrix_mask_str(m))

    n = m.nrows()
    deps_1_1 = {}
    for i, row in enumerate(m.rows()):
        deps_1_1[i] = set(j for j,ex in enumerate(row) if not bool(ex==0))

    deps_1_all = {}
    def find_deps_1_all(i):
        if deps_1_all.has_key(i):
            return deps_1_all[i]
        deps_1_all[i] = deps_1_1[i]
        for j in deps_1_1[i]:
            if i == j:
                continue
            find_deps_1_all(j)
            deps_1_all[i] = deps_1_all[i].union(deps_1_all[j])
    [find_deps_1_all(j) for j in xrange(n)]

    deps_coup = dict((i, set([])) for i in xrange(n))
    for i in xrange(n):
        for j in xrange(n):
            if (i in deps_1_all[j]) and (j in deps_1_all[i]):
                deps_coup[i].update([i,j])
                deps_coup[j].update([i,j])

    shuffle = []
    error = False
    while not error:
        if not deps_coup:
            break
        error = True
        for i in xrange(n):
            if not deps_coup.has_key(i):
                continue
            b = deps_coup[i].copy()
            if b != deps_1_all[i]:
                continue
            if b == set([]):
                b = set([i])
            shuffle += [b]
            for j in deps_1_all[i].copy():
                if deps_1_all.has_key(j):
                    del deps_1_all[j]
                    del deps_coup[j]
            for j in xrange(n):
                if not deps_coup.has_key(j):
                    continue
                deps_coup[j].difference_update(b)
                deps_1_all[j].difference_update(b)
            if deps_coup.has_key(i):
                del deps_coup[i]
            if deps_1_all.has_key(i):
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
    mt = transform(m, None, t)
    if logger.isEnabledFor(logging.INFO):
        logger.info("Matrix mask after transformation:\n%s\n" % matrix_mask_str(mt))

    return mt, t, blocks

def canonical_form(m, x, eps):
    m, t1, b = block_triangular_form(m)
    m, t2 = normalize_by_blocks(m, b, x, eps)
    m, t3 = fuchsify_by_blocks(m, b, x, eps)
    m, t4 = factorize(m, x, eps)
    t = t1*t2*t3*t4
    return m, t

def matrix_mask(m):
    n = m.nrows()
    return matrix(SR, n, n, [int(not bool(ex==0)) for ex in m.list()])

def matrix_mask_str(m):
    s = ''
    for row in matrix_mask(m).rows():
        s += ' '.join([str(ex) for ex in row]).replace('0', '.').replace('1','x') + "\n"
    return s

#==================================================================================================
# Step I: Fuchsify
#==================================================================================================

def is_fuchsian(m, x, cas="maxima"):
    for i, expr in enumerate(m.list()):
        if expr.is_zero():
            continue
        points = singularities_expr(expr, x, cas=cas)
        for x0, p in points.iteritems():
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
    assert M.is_square()
    rng = Random(seed)
    poincare_map = singularities(M, x)
    def iter_reductions(p1, U):
        for p2, prank2 in poincare_map.iteritems():
            if bool(p2 == p1): continue
            while prank2 >= 0:
                B0 = matrix_c0(M, x, p2, prank2)
                if not B0.is_zero(): break
                poincare_map[p2] = prank2 = prank2 - 1
            if prank2 < 0: continue
            v = find_dual_basis_spanning_left_invariant_subspace(B0, U, rng)
            if v is None: continue
            P = (U*v).simplify_rational()
            M2 = balance_transform(M, P, p1, p2, x).simplify_rational()
            yield p2, P, M2
    combinedT = identity_matrix(M.base_ring(), M.nrows())
    reduction_points = [pt for pt,p in poincare_map.iteritems() if p >= 1]
    reduction_points.sort()
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
                        "%s\nfurther reduction is pointless:\n%s",
                        M, e)
                raise FuchsiaError("matrix cannot be reduced to Fuchsian form")
            try:
                point2, P, M = min(iter_reductions(point, U), \
                        key=lambda (point2, P, M2): matrix_complexity(M2))
            except ValueError as e:
                point2 = any_integer(rng, M.base_ring(), poincare_map)
                P = (U*V).simplify_rational()
                M = balance_transform(M, P, point, point2, x)
                M = M.simplify_rational()
                logger.info(
                        "Will introduce an apparent singularity at %s.",
                        point2)
            logger.debug(
                    "Applying balance between %s and %s with projector:\n%s",
                    point, point2, P)
            combinedT = combinedT * balance(P, point, point2, x)
            if point2 not in poincare_map:
                poincare_map[point2] = 1
        poincare_map[point] = prank - 1
    combinedT = combinedT.simplify_rational()
    return M, combinedT

def fuchsify_by_blocks(m, b, x, eps, cas="maxima"):
    n = m.nrows()
    m0, t = m, identity_matrix(SR, n)
    for i, (ki, ni) in enumerate(b):
        for j, (kj, nj) in enumerate(reversed(b[:i])):
            logger.info("processing blocks (%d, %d) (%d, %d)\n" % (kj, nj, ki, ni))
            pts = singularities(m.submatrix(ki, kj, ni, nj), x, cas=cas)
            logger.info("singular points = %s" % pts)
            while any(pts.values()):
                bj = m.submatrix(ki, kj, ni, nj)
                if bj.is_zero():
                    break
                for x0, p in pts.iteritems():
                    if p < 1:
                        continue
                    a0 = matrix_residue(m.submatrix(ki, ki, ni, ni)/eps, x, x0)
                    b0 = matrix_c0(bj, x, x0, p, cas=cas)
                    c0 = matrix_residue(m.submatrix(kj, kj, nj, nj)/eps, x, x0)
        
                    d_vars = [gensym() for i in xrange(ni*nj)]
                    d = matrix(SR, ni, nj, d_vars)
                    eq = d + eps/p*(a0*d - d*c0) + b0/p
                    sol = solve(eq.list(), *d_vars, solution_dict=True)
                    d = d.subs(sol[0])

                    t0 = identity_matrix(SR, n)
                    t0[ki:ki+ni, kj:kj+nj] = \
                            d/(x-x0)**p if not (x0 == oo) else -d*(x**p)
                    m = transform(m, x, t0).simplify_rational()

                    t = (t*t0).simplify_rational()
                    pts[x0] -= 1
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
        M = M.simplify_rational()
        combinedT = combinedT * balance(P, v, v2, x)
    combinedT = combinedT.simplify_rational()
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
    assert all(A0J_cs[i] >= A0J_cs[i+1] for i in xrange(len(A0J_cs) - 1))
    ncells = len(A0J_cs)
    A0J_css = [sum(A0J_cs[:i]) for i in xrange(ncells + 1)]
    nsimplecells = sum(1 if s == 1 else 0 for s in A0J_cs)
    u0 = [U[:,A0J_css[i]] for i in xrange(ncells)]
    v0t = [invU[A0J_css[i+1]-1,:] for i in xrange(ncells)]
    L0 = matrix([
        [(v0t[k]*(A1)*u0[l])[0,0] for l in xrange(ncells)]
        for k in xrange(ncells)
    ])
    L0 = L0.simplify_rational()
    L1 = matrix([
        [(v0t[k]*u0[l])[0,0] for l in xrange(ncells)]
        for k in xrange(ncells)
    ])
    assert (L1 - diagonal_matrix(
            [0]*(ncells-nsimplecells) + [1]*nsimplecells)).is_zero()
    #zero_rows = [i for i in xrange(A0.nrows()) if A0J[i,:].is_zero()]
    #zero_cols = [j for j in xrange(A0.nrows()) if A0J[:,j].is_zero()]
    #assert len(zero_rows) == len(zero_cols) == ncells
    #_L0 = (invU*A1*U)[zero_rows, zero_cols]
    #assert (L0 - _L0).is_zero()
    lam = SR.symbol()
    if not (L0 - lam*L1).determinant().is_zero():
        raise FuchsiaError("matrix is Moser-irreducible")
    S, D = alg1(L0, A0J_cs)
    I_E = identity_matrix(D.base_ring(), A0.nrows())
    for i in xrange(ncells):
        for j in xrange(ncells):
            if not D[i,j].is_zero():
                ni = A0J_css[i]
                nj = A0J_css[j]
                for k in xrange(min(A0J_cs[i], A0J_cs[j])):
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
        for i in xrange(0, N):
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
        for j in xrange(fi):
            D0[j, fi] = -c[j, 0]
            invD0[j, fi] = -c[j, 0] \
                    if jordan_cellsizes[j] == jordan_cellsizes[fi] else 0
        L0 = (I_N - invD0)*L0*(I_N + D0)
        D = D + D0 + D*D0
        if fi < K:
            break
        S.add(fi)
    # Check Alg.1 promise
    for j in xrange(N):
        for k in xrange(N):
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
    for x0, p in points.iteritems():
        M0 = matrix_residue(M, x, x0)
        for ev in M0.eigenvalues():
            ev = limit_fixed(ev, eps, 0)
            if not (Rational((-1, 2)) <= ev and ev < Rational((1, 2))):
                return False
    return True

def is_normalized_by_blocks(m, b, x, eps):
    """Return True if diagonal blocks of `m` are normalized, that is all the eigenvalues of theirs
    residues in `x` lie in the range [-1/2, 1/2) in eps->0 limit; return False otherwise. Diagonal
    blocks are defined by the list `b` which corresponds to the equivalent value returned by the
    `block_triangular_form` routine.

    Examples:
    >>> x, e = var("x epsilon")
    >>> is_normalized_by_blocks(matrix([[(1+e)/3/x, 0], [0, e/x]]), [(0,1),(1,1)], x, e)
    True
    """
    for ki, ni in b:
        mi = m.submatrix(ki, ki, ni, ni).simplify_rational()
        if not is_normalized(mi, x, eps):
            return False
    return True

def normalize(m, x, eps, seed=0):
    """Given a Fuchsian system of differential equations of the
    form dF/dx=m*F, find a transformation that will shift all
    the eigenvalues of m's residues into [-1/2, 1/2) range (in
    the limit eps->0). Return the transformed matrix m and the
    transformation. Raise FuchsiaError if such transformation
    is not found.
    """
    logger.info("[normalize] START\n")
    T = identity_matrix(m.base_ring(), m.nrows())
    state = {"random": Random(seed)}
    i = 0
    while True:
        i += 1
        m = partial_fraction(m, x)
        logger.info("[normalize] STEP #%s\n" % i)
        balances = find_balances(m, x, eps, state)
        b = select_balance(balances, eps, state)
        if b is None:
            if state["is_normalized"]:
                break
            raise FuchsiaError("can not balance matrix")
        logger.info("Use the balance:\n    %s\n" % b)

        cond, x1, x2, a0_eval, b0_eval, a0_evec, b0_evec, scale = b
        if cond == 1:
            P = cross_product(a0_evec, b0_evec) / scale
            m = balance_transform(m, P, x1, x2, x)
            T0 = balance(P, x1, x2, x)
        else:
            P = cross_product(b0_evec, a0_evec) / scale
            m = balance_transform(m, P, x2, x1, x)
            T0 = balance(P, x2, x1, x)

        T = (T*T0).simplify_rational()
    logger.info("[normalize] DONE\n")
    return m.simplify_rational(), T

def normalize_by_blocks(m, b, x, eps):
    """Given a lower block-triangular system of differential equations of the form dF/dx=m*F,
    find a transformation that will shift all eigenvalues of all residues of all its diagonal
    blocks into the range [-1/2, 1/2) in eps->0 limit. Return the transformed matrix m and the
    transformation; raise FuchsiaError if such transformation is not found. Diagonal blocks
     are defined by the list `b` which corresponds to the equivalent value returned by the
    `block_triangular_form` routine.
    """
    n = m.nrows()
    t = identity_matrix(SR, n)
    for i, (ki, ni) in enumerate(b):
        mi = m.submatrix(ki, ki, ni, ni).simplify_rational()
        ti = identity_matrix(SR, ni)
        if is_verbose():
            msg = ("Reducing block #%d (%d,%d):\n%s\n" % (i,ki,ni,mi)).replace("\n", "\n  ")
            logger.info(msg)

        mi_fuchs, ti_fuchs = fuchsify(mi, x, eps)
        ti = ti*ti_fuchs

        mi_norm, ti_norm = normalize(mi_fuchs, x, eps)
        ti = ti*ti_norm

        mi_eps, ti_eps = factorize(mi_norm, x, eps)
        ti = ti*ti_eps

        t[ki:ki+ni, ki:ki+ni] = ti
    mt = transform(m, x, t).simplify_rational()
    return mt, t

def find_balances(m, x, eps, state={}):
    if not state.has_key("pairs"):
        pairs = list(combinations(singularities(m, x).keys(), 2))
        state["pairs"] = []
        for x1,x2 in pairs:
            state["pairs"] += [(x1,x2),(x2,x1)]
            # the order of pairs matters, hence we do not use `permutations` routine
    pairs = state["pairs"]
    pairs.sort(key = lambda pair: pair in state.get("processed_pairs", []))
    state["is_normalized"] = True
    state.setdefault("processed_pairs", [])

    for x1, x2 in pairs:
        logger.info("Looking for the balance between x = %s and x = %s" % (x1,x2))
        a0, b0 = matrix_residue(m, x, x1, cas="maxima"), matrix_residue(m, x, x2, cas="maxima")

        a0_evr, b0_evl = eigenvectors_right(a0), eigenvectors_left(b0)
        for ev, evec, emult in a0_evr:
            if limit_fixed(ev, eps, 0) != 0:
                state["is_normalized"] = False
                break

        if logger.isEnabledFor(logging.INFO):
            msg = "  Eigenvalues:\n"
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
            logger.info(msg)

        balances_1 = find_balances_by_cond(a0_evr, b0_evl, lambda a0_eval, b0_eval: limit_fixed(a0_eval, eps, 0) < -0.5)
        for balance in balances_1:
            balance = [1, x1, x2] + balance
            yield balance

        a0_evl, b0_evr = eigenvectors_left(a0), eigenvectors_right(b0)
        balances_2 = find_balances_by_cond(a0_evl, b0_evr, lambda a0_eval, b0_eval: limit_fixed(a0_eval, eps, 0) >= 0.5)
        for balance in balances_2:
            balance = [2, x1, x2] + balance
            yield balance

        if len(state["processed_pairs"]) == len(pairs):
            state["processed_pairs"] = []
        if (x1,x2) not in state["processed_pairs"]:
            state["processed_pairs"].append((x1,x2))

def find_balances_by_cond(a0_ev, b0_ev, cond):
    res = []
    for a0_eval, a0_evecs, a0_evmult in a0_ev:
        for b0_eval, b0_evecs, b0_evmult in b0_ev:
            if not cond(a0_eval, b0_eval):
                logger.debug("Balance rejected:\n    a0_eval = %s\n    b0_eval = %s" % (a0_eval, b0_eval))
                continue
            for a0_evec in a0_evecs:
                for b0_evec in b0_evecs:
                    scale = dot_product(a0_evec, b0_evec).simplify_rational()
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

    x0 = state.get("x0")
    if x0 is None:
        for b in bs:
            cond, x1, x2, ev1, ev2 = b[:5]
            if cond == 1:
                x0 = x2
                break
            if cond == 2:
                x0 = x1
                break
        logger.info("Select x0 = %s" % x0)
        state["x0"] = x0

    balances_x0 = [b for b in bs if (b[0] == 1 and b[2] == x0) or (b[0] == 2 and b[1] == x0)]
    b = state["random"].choice(balances_x0) if balances_x0 else None
    return b

def eigenvectors_left(m):
    if m._cache is None:
        m._cache = {}
    if not m._cache.has_key("eigenvectors_left"):
        res = simplify(m.eigenvectors_left())
        m._cache["eigenvectors_left"] = res
    return m._cache["eigenvectors_left"]

def eigenvectors_right(m):
    if m._cache is None:
        m._cache = {}
    if not m._cache.has_key("eigenvectors_right"):
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

def f_solve(eqs, var, cas="maxima"):
    if cas == "maxima":
        return solve(eqs, var, solution_dict=True)
    elif cas == "maple":
        s = maple.solve(eqs, var)
        solution = s.parent().get(s._name).split(',')
        if solution == ['']:
            return []
        res = []
        for s in solution:
            try:
                expr = expand(simplify(_parser.parse(s)))
                res.append(expr)
            except SyntaxError as error:
                print "ERROR:  \n%s\n  %s\n" % (s, error)
                continue
        return res
    raise FuchsiError("Fuchsia does not have support for 'solve' routine of %s" % cas)

def factorize(M, x, epsilon, seed=0):
    """Given a normalized Fuchsian system of differential equations:
        dF/dx = M(x,epsilon)*F,
    try to find a transformation that will factor out an epsilon
    from M. Return a transformed M (proportional to epsilon)
    and T. Raise FuchsiaError if epsilon can not be factored.
    """
    n = M.nrows()
    M = M.simplify_rational()
    if epsilon not in (M/epsilon).variables():
        return M, identity_matrix(SR, n)
    rng = Random(seed)
    mu = gensym()
    T_symbols = [gensym() for i in xrange(n*n)]
    T = matrix(SR, n, n, T_symbols)
    eqs = []
    for point, prank in singularities(M, x).iteritems():
        assert prank == 0
        R = matrix_c0(M, x, point, 0)
        eq = (R/epsilon)*T-T*(R.subs({epsilon: mu})/mu)
        eqs.extend(eq.list())
    for solution in solve(eqs, T_symbols, solution_dict=True):
        S = T.subs(solution)
        if not S.is_invertible():
            continue
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
                sT = sT.simplify_rational()
                M = transform(M, x, sT).simplify_rational()
            except (ZeroDivisionError, ValueError):
                rndrange += 1 + rndrange//4
                continue
            break
        # We're leaking a bunch of temprary variables here,
        # which accumulate in SR.variables, but do I care?
        return M, sT
        # No, I don't.
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
    for point, prank in singularities(M, x).iteritems():
        R = matrix_c0(M, x, point, prank)
        J, T = R.jordan_form(transformation=True)
        MM = transform(M, x, T).simplify_rational()
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
        for factor, n in factors.iteritems():
            if not filter(factor): continue
            if factor in factor2exp2count:
                factor2exp2count[factor][n] += 1
            else:
                if i > 0: factor2exp2count[factor][0] = i
                factor2exp2count[factor][n] = 1
        for factor, exps in factor2exp2count.iteritems():
            if factor not in factors:
                exps[0] += 1
    result = SR(1)
    for factor, exp2count in factor2exp2count.iteritems():
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
    for i in xrange(n):
        factor = common_factor(
            [M[i,k] for k in xrange(n) if i != k and not M[i,k].is_zero()] +
            [1/M[k,i] for k in xrange(n) if k != i and not M[k,i].is_zero()],
            lambda e: x not in e.variables())
        if factor != 1:
            dT = identity_matrix(M.base_ring(), n)
            dT[i,i] = T[i,i] = factor
            M = transform(M, x, dT)
    logger.info(
            "stripping common factors with this transform:\n"
            "diagonal_matrix([\n  %s\n])",
            ",\n  ".join(str(e) for e in T.diagonal()))
    return M.simplify_rational(), T

#==============================================================================
# Import/Export routines
#==============================================================================

def import_matrix_from_file(filename, fmt="mtx"):
    """Read and return a matrix stored in a given format from a named file.
    """
    with open(filename, 'r') as f:
        if fmt == "mtx":
            return import_matrix_matrixmarket(f)
        elif fmt == "m":
            return import_matrix_mathematica(f)
        else:
            raise FuchsiaError("Unknown matrix format '%s'" % fmt)

_parser = Parser(make_int=ZZ, make_float=RR, make_var=SR.var)

def import_matrix_mathematica(f):
    """Read and return a matrix, stored in the Mathematica format, from a file-like object.
    """
    data =  f.read()
    data = data.translate(None, ' \n')
    nr = data.count('},{')+1
    data = data.translate(None, '{}')
    data = data.replace(',', '\n')
    nc = int((data.count('\n')+1)/nr)

    import StringIO
    sio = StringIO.StringIO("%d %d\n" % (nc, nr) + data)
    return import_matrix_matrixmarket(sio).T

def import_matrix_matrixmarket(f):
    """Read and return a matrix, stored in the MatrixMarket array format, from a file-like object.
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
        msg = ("Can not determine size of the matrix in '%s'\n"
               "Make sure the file is in the MatrixMarket array format")
        raise FuchsiaError(msg % f.name)
    data = []
    for s in f.readlines():
        lineno += 1
        if s.startswith('%'):
            continue
        try:
            expr = _parser.parse(s)
        except SyntaxError:
            msg = ("Malformad expression at line %d of '%s':\n%s\n"
                   "Make sure the file is in the MatrixMarket array format")
            raise FuchsiaError(msg % (lineno, f.name, s))
        data.append(expr)
    m = matrix(ncols, nrows, data).T
    return m

def export_matrix_to_file(filename, m, fmt="mtx"):
    """Write a matrix in a given format to a named file.
    """
    with open(filename, 'w') as f:
        if fmt == "mtx":
            export_matrix_matrixmarket(f, m)
        elif fmt == "m":
            export_matrix_mathematica(f, m)
        else:
            raise FuchsiaError("Unknown matrix format '%s'" % fmt)

def export_matrix_mathematica(f, m):
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
    """Write a matrix to a file-like object using a given format.
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
        print __doc__
        exit(0)

    try:
        print('\033[35;1mFuchsia v%s (commit: %s)\033[0m' % (__version__, __commit__))
        print "Authors: %s\n" % __author__
        mpath = tpath = profpath = fmt = None
        M = T = x_fy = None
        x, epsilon = SR.var("x eps")
        fmt = "mtx"
        logger.setLevel(logging.INFO)
        kwargs, args = getopt.gnu_getopt(sys.argv[1:], "hvl:f:P:x:y:e:m:t:")
        for key, value in kwargs:
            if key == "-h": usage()
            if key == "-f": fmt = value
            if key == "-l":
                fh = logging.FileHandler(value, "w")
                fh.setFormatter(logging.Formatter(logger_format))
                logger.addHandler(fh)
            if key == "-v": logger.setLevel(logging.DEBUG)
            if key == "-P": profpath = value
            if key == "-x": x = _parser.parse(value)
            if key == "-y": x_fy = _parser.parse(value)
            if key == "-e": epsilon = SR.var(value)
            if key == "-m": mpath = value
            if key == "-t": tpath = value
        with profile(profpath):
            if len(args) == 2 and args[0] == 'fuchsify':
                M = import_matrix_from_file(args[1], fmt=fmt)
                M, t1 = simplify_by_factorization(M, x)
                M, t2 = fuchsify(M, x)
                T = t1*t2
            elif len(args) == 2 and args[0] == 'normalize':
                M = import_matrix_from_file(args[1], fmt=fmt)
                M, t1 = simplify_by_factorization(M, x)
                i = 0
                M, t2 = normalize(M, x, epsilon)
            elif len(args) == 2 and args[0] == 'factorize':
                M = import_matrix_from_file(args[1], fmt=fmt)
                M, t1 = simplify_by_factorization(M, x)
                M, t2 = factorize(M, x, epsilon)
                T = t1*t2
            elif len(args) == 2 and args[0] == 'sort':
                M = import_matrix_from_file(args[1], fmt=fmt)
                M, T, B = block_triangular_form(M)
            elif len(args) == 2 and args[0] == 'reduce':
                m = import_matrix_from_file(args[1], fmt=fmt)
                M, T = canonical_form(m, x, epsilon)
            elif len(args) >= 2 and args[0] == 'transform':
                M = import_matrix_from_file(args[1], fmt=fmt)
                if x_fy is None:
                    t = import_matrix_from_file(args[2], fmt=fmt)
                    M = transform(M, x, t)
                else:
                    x, fy = x_fy.left(), x_fy.right()
                    y = None
                    [y for y in fy.variables() if y != epsilon]
                    if y is None:
                        raise FuchsiaError("No new variable found in the expression given by option -y '%s'" % x_fy)
                    assert y.is_symbol()
                    M = change_variable(M, x, y, fy)
                    x = y
            elif len(args) == 0:
                usage()
            else:
                raise getopt.GetoptError("unknown command: %s" % (args))
        if M is not None:
            M = partial_fraction(M, x)
            if mpath is not None: 
                export_matrix_to_file(mpath, M, fmt=fmt)
            else:
                print M
        if tpath is not None and T is not None:
            T = partial_fraction(T, x)
            export_matrix_to_file(tpath, T, fmt=fmt)
    except getopt.GetoptError as error:
        logger.error("%s", error)
        exit(1)

if __name__ == '__main__':
    try:
        main()
    except FuchsiaError as e:
        msg, tab = str(e), '  '
        msg = tab+msg.replace('\n', '\n'+tab)
        logger.error(msg)
        exit(1)
    except Exception as e:
        if is_verbose():
            raise
