"""
    References:
        [1] Roman Lee, arXiv:1411.0911
"""
from   collections import defaultdict
from   itertools import permutations
import logging
import os.path as path
from   random import Random

from   sage.all import *
from   sage.misc.parser import Parser

logging.basicConfig(
    format='\033[32m%(levelname)s [%(asctime)s]\033[0m\n%(message)s',
#    format='\033[32m%(levelname)s [%(asctime)s]\033[0m %(pathname)s:%(funcName)s, line %(lineno)s %(message)s\n',
#    format='%(levelname)s [%(asctime)s] %(message)s %(pathname)s:%(funcName)s, line %(lineno)s',
    datefmt='%Y-%m-%d %I:%M:%S',
)
logger = logging.getLogger('fuchsia')
DEBUG = logger.getEffectiveLevel() <= logging.DEBUG
INFO  = logger.getEffectiveLevel() <= logging.INFO

__author__ = "Oleksandr Gituliar, Vitaly Magerya"
__author_email__ = "oleksandr.gituliar@ifj.edu.pl"
__version__ = open(path.join(path.dirname(__file__),"VERSION")).readline().rstrip('\n')
try:
    __commit__ = open(path.join(path.dirname(__file__),"COMMIT")).readline().rstrip('\n')
except IOError:
    __commit__ = "unknown"

def partial_fraction(M, var):
    return M.apply_map(lambda ex: ex.partial_fraction(var))

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
    assert x1 != x2
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
    assert x1 != x2
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

def limit_fixed(expr, x, lim):
    """Return a limit of expr when x->lim.

    The standard 'limit()' function of SageMath does not allow
    you to specify the variable, only it's name as a keyword
    argument. If you have a variable, and not it's name, use
    this function instead.
    """
    l = maxima_calculus.sr_limit(expr, x, lim)
    return expr.parent()(l)

def poincare_rank_at_oo(M, x):
    """Return Poincare rank of matrix M at x=Infinity.

    Examples:
    >>> x = var("x")
    >>> poincare_rank_at_oo(matrix([[x, 1, 1/x]]), x)
    2
    >>> poincare_rank_at_oo(matrix([[1]]), x)
    1
    >>> poincare_rank_at_oo(matrix([[1/x]]), x)
    0
    """
    p = -1
    for expr in (M.subs({x: 1/x})/x**2).list():
        if x not in expr.variables():
            continue
        solutions, ks = solve(Integer(1)/expr, x,
                solution_dict=True, multiplicities=True)
        for sol, k in zip(solutions, ks):
            if sol[x] == 0:
                p = max(p, k - 1)
    return p

def singularities(M, x):
    """Find values of x around which rational matrix M has
    a singularity; return a dictionary with {val: p} entries,
    where p is the Poincare rank of M at x=val.

    Example:
    >>> x, y = var("x y")
    >>> M = matrix([[1/x, 0], [1/(x+y)**3, 1/(x+y)]])
    >>> s = singularities(M, x)
    >>> from pprint import pprint; pprint(s)
    {0: 0, -y: 2, +Infinity: 0}
    """
    M = M.simplify_rational()
    result = {}
    for expr in M.list():
        if x not in expr.variables():
            # If expression is constant in x, no singularity is present.
            # In particular, this prevents solve(1/0) calls.
            continue
        solutions, ks = solve(Integer(1)/expr, x,
                solution_dict=True, multiplicities=True)
        for sol, k in zip(solutions, ks):
            val = sol[x]
            if val in result:
                result[val] = max(result[val], k - 1)
            else:
                result[val] = k - 1
    p = poincare_rank_at_oo(M, x)
    if p >= 0:
        result[oo] = p
    return result

def matrix_taylor0(M, x, point, exp):
    """Return the 0-th coefficient of Taylor expansion of
    a matrix M around a finite point x=point, assuming that
    M(x->point)~1/(x-point)**exp.

    Example:
    >>> x = var('x')
    >>> matrix_taylor0(matrix([[x/(x-1), 1/x, x, 1]]), x, 1, 1)
    [1 0 0 0]
    """
    return matrix([
        [taylor(e, x, 0, 0) for e in row]
        for row in M.subs({x: x+point})*x**exp
    ])

def matrix_taylor1(M, x, point, exp):
    """Return the 1-th coefficient of Taylor expansion of
    a matrix M around a finite point x=point, assuming that
    M(x->point)~1/(x-point)**exp.

    Example:
    >>> x = var('x')
    >>> matrix_taylor1(matrix([[x/(x-1), 1/x, x, 1]]), x, 0, 1)
    [0 0 0 1]
    """
    return matrix([
        [taylor(e, x, 0, 1).coefficient(x) for e in row]
        for row in M.subs({x: x+point})*x**exp
    ])

def matrix_c0(M, x, point, p):
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
        return matrix_taylor0(M.subs({x: 1/x}), x, 0, p-1)
    else:
        return matrix_taylor0(M, x, point, p+1)

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
    m0 = matrix_c0(M, x, x0, 0)
    if x0 == oo:
        return -m0
    else:
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
        raise ValueError, "matrix equation has no solutions"
    return C

def solve_left_fixed(A, B):
    """As of SageMath 6.10, 'Matrix.solve_left' method uses a broken
    check for solution correctness; this function corrects that check.
    """
    return solve_right_fixed(A.transpose(), B.transpose()).transpose()

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

def alg1x(A0, A1, x):
    if not matrix_is_nilpotent(A0):
        raise ValueError("matrix is irreducible (non-nilpotent residue)")
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
        raise ValueError("matrix is Moser-irreducible")
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

def any_integer(rng, ring, excluded):
    r = 2
    while True:
        p = ring(rng.randint(-r, r))
        if p not in excluded:
            return p
        r *= 2

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

def fuchsify(M, x, seed=0):
    """Given a system of differential equations of the form dF/dx=M*F,
    try to find a transformation T, which will reduce M to Fuchsian
    form. Return the transformed M and T.

    Note that such transformations are not unique; you can obtain
    different ones by supplying different seeds.
    """
    assert M.is_square()
    rng = Random(seed)
    poincare_map = singularities(M, x)
    def iter_reductions(p1, U):
        for p2, prank2 in poincare_map.iteritems():
            if p2 == p1: continue
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
            U, V = alg1x(A0, A1, x)
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

def normalize(m, x, eps, seed=0):
    m = partial_fraction(m, x)
    T = identity_matrix(m.base_ring(), m.nrows())
    select_balance_state = {"random": Random(seed)}
    while not is_normalized(m, x, eps):
        if INFO:
            points = singularities(m, x).keys()
            msg = "Eigenvalues:\n"
            for x0 in points:
                m0 = matrix_residue(m, x, x0)
                msg += "    x = %s:\n" % x0
                msg += "        %s\n" % str(m0.eigenvalues()).replace("\n"," ")
                msg += "        l: %s\n" % str(m0.eigenvectors_left()).replace("\n"," ")
                msg += "        r: %s\n" % str(m0.eigenvectors_right()).replace("\n"," ")
            logger.info(msg)

        balances = find_balances(m, x, eps)
        b = select_balance(balances, eps, select_balance_state)
        if b is None:
            logger.info("Can not balance matrix")
            raise ValueError
        logger.info("Use balance:\n    %s" % b)

        cond, x1, x2, a0_eval, b0_eval, a0_evec, b0_evec, scale = b
        if cond == 1:
            P = partial_fraction(cross_product(a0_evec, b0_evec) / scale, eps)
            T0 = balance(P, x1, x2, x)
        else:
            P = partial_fraction(cross_product(b0_evec, a0_evec) / scale, eps)
            T0 = balance(P, x2, x1, x)
        T0 = partial_fraction(T0, x)
        logger.debug("P  =\n%s" % (P,))
        logger.debug("T0 =\n%s" % (T0,))

        m = transform(m, x, T0)
        m = partial_fraction(m, x)
        if INFO:
            logger.info("New matrix:\n    %s" % '\n    '.join([str(ex) for ex in m.list()]))
        T = partial_fraction(T*T0, x)

    logger.info("Transformation matrix:\n    %s" % '\n    '.join([str(ex) for ex in T.list()]))
    return m, T

def find_balances(m, x, eps):
    points = singularities(m, x).keys()
    residues = [matrix_residue(m, x, pt) for pt in points]
    left_evectors = [r.eigenvectors_left() for r in residues]
    right_evectors = [r.eigenvectors_right() for r in residues]
    balances = []
    for i1, i2 in permutations(range(len(points)), 2):
        x1, x2 = points[i1], points[i2]
        a0, b0 = residues[i1], residues[i2]
        a0_evr, b0_evr = right_evectors[i1], right_evectors[i2]
        a0_evl, b0_evl = left_evectors[i1], left_evectors[i2]
        logger.debug("Looking for the balance\n    x1 = %s\n    x2 = %s" % (x1, x2))
        if DEBUG:
            logger.debug("Res(m, x = %s)\n%s" % (x1, a0))
            print("\nRes(m, x = %s)\n%s" % (x2, b0))
            logger.debug("Eigenvalues:")
            for xi, mi in [(x1, a0), (x2, b0)]:
                print("    x = %s:" % xi)
                print("        %s" % str(mi.eigenvalues()).replace("\n"," "))
                print("        l: %s" % str(mi.eigenvectors_left()).replace("\n"," "))
                print("        r: %s" % str(mi.eigenvectors_right()).replace("\n"," "))


        logging.debug("Use condition #1")
        balances_1 = find_balances_by_cond(a0_evr, b0_evl, lambda a0_eval, b0_eval: limit_fixed(a0_eval, eps, 0) < -0.5)
        balances_1 = [[1, x1, x2] + balance for balance in balances_1]
        balances += balances_1
    
        logging.debug("Use condition #2")
        balances_2 = find_balances_by_cond(a0_evl, b0_evr, lambda a0_eval, b0_eval: limit_fixed(a0_eval, eps, 0) >= 0.5)
        balances_2 = [[2, x1, x2] + balance for balance in balances_2]
        balances += balances_2

    if INFO:
        logger.info("Found balances (#cond, x1, x2, a0_eval, b0_eval, a0_evec, b0_evec, scale):")
        for balance in balances:
            print "    %s" % (balance,)
    return balances

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

def select_balance(balances, eps, state):
    for b in balances:
        cond, x1, x2, a0_eval, b0_eval, a0_evec, b0_evec, scale = b
        if (cond == 1) and limit_fixed(a0_eval, eps, 0) < -0.5 and \
                limit_fixed(b0_eval, eps, 0) >= 0.5:
            return b
        elif (cond == 2) and limit_fixed(a0_eval, eps, 0) >= 0.5 and \
                limit_fixed(b0_eval, eps, 0) < -0.5:
            return b

    x0 = state.get("x0")
    if x0 is None:
        for b in balances:
            cond, x1, x2, ev1, ev2 = b[:5]
            if cond == 1:
                x0 = x2
                break
            if cond == 2:
                x0 = x1
                break
        logger.info("Select x0 = %s" % x0)
        state["x0"] = x0

    balances_x0 = [b for b in balances if (b[0] == 1 and b[2] == x0) or (b[0] == 2 and b[1] == x0)]
    b = state["random"].choice(balances_x0) if balances_x0 else None
    return b

def gensym():
    sym = SR.symbol()
    SR.symbols[str(sym)] = sym
    return sym

def factor_epsilon(M, x, epsilon, seed=0):
    """Given a normalized Fuchsian system of differential equations:
        dF/dx = M(x,epsilon)*F,
    try to find a transformation that will factor out an epsilon
    from M. Return a transformed M (proportional to epsilon) and T.
    """
    n = M.nrows()
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
        while True:
            try:
                sT = S.subs([
                    e==rng.randint(-99, 99)
                    for e in S.variables() if e != epsilon
                ])
                sT = sT.simplify_rational()
                M = transform(M, x, sT).simplify_rational()
            except (ZeroDivisionError, ValueError):
                continue
            break
        # We're leaking a bunch of temprary variables here,
        # which accumulate in SR.variables, but do I care?
        return M, sT
        # No, I don't.
    raise ValueError("can not factor epsilon")

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
    return M.simplify_rational(), T

def is_normalized(M, x, eps):
    """Return True if (a Fuchsian) matrix M is normalized, that
    is all the eigenvalues of it's residues lie in [-1/2, 1/2)
    range (in limit eps->0). Return False otherwise.

    Examples:
    >>> x, e = var("x epsilon")
    >>> is_normalized(matrix([[(1+e)/3/x, 0], [0, e/x]]), x, e)
    True
    """
    points = singularities(M, x)
    for x0, p in points.iteritems():
        # m should also be Fuchsian
        if p > 0:
            return False
        M0 = matrix_residue(M, x, x0)
        for ev in M0.eigenvalues():
            ev = limit_fixed(ev, eps, 0)
            if not (Rational((-1, 2)) <= ev and ev < Rational((1, 2))):
                return False
    return True

# matrix.py

def export_matrix(f, m):
    """Write a matrix to a file-like object using MatrixMarket
    array format.
    """
    f.write("%%MatrixMarket matrix array Maple[symbolic] general\n")
    f.write("%d %d\n" % (m.nrows(), m.ncols()))
    for col in m.columns():
        for mij in col:
            f.write(str(mij).replace(' ', ''))
            f.write('\n')

def export_matrix_to_file(filename, m):
    """Write a matrix to a named file using MatrixMarket array format.
    """
    with open(filename, 'w') as f:
        export_matrix(f, m)

_parser = Parser(make_int=ZZ, make_float=RR, make_var=SR.var)

def import_matrix(f):
    """Read and return a matrix, stored in the MatrixMarket
    array format, from a file-like object.
    """
    while True:
        s = f.readline()
        if not s.startswith('%'):
            break
    nrows, ncols = map(int, s.split())
    data = [_parser.parse(s) for s in f.readlines() if not s.startswith('%')]
    m = matrix(ncols, nrows, data).T
    return m

def import_matrix_from_file(filename):
    """Read and return a matrix stored in MatrixMarket array
    format from a named file.
    """
    with open(filename, 'r') as f:
        return import_matrix(f)

def cross_product(v1, v2):
    m1, m2 = matrix(v1), matrix(v2)
    return m1.transpose() * m2

def dot_product(v1, v2):
    m1, m2 = matrix(v1), matrix(v2)
    sp = m1 * m2.transpose()
    return sp[0,0]
