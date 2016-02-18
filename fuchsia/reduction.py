"""
    References:
        [1] Roman Lee, arXiv:1411.0911
"""
from   itertools import combinations
import logging

from   sage.all import *
from   random import Random

from   fuchsia.matrix import (cross_product, dot_product)

logger = logging.getLogger('fuchsia')

def transform(M, x, T):
    """Given a system of differential equations dF/dx = M*F,
    and a transformation of base functions F = T*F', compute
    and return M', such that dF'/dx = M'*F'.

    Note: M' = inverse(T)*(M*T - dT/dx)
    """
    return T.inverse()*(M*T - derivative(T, x))

def balance(P, x1, x2, x):
    assert P.is_square()
    #assert (P*P - P).is_zero()
    assert x1 != x2
    coP = identity_matrix(P.nrows()) - P
    if x1 == oo:
        return coP - (x - x2)*P
    elif x2 == oo:
        return coP - 1/(x - x1)*P
    else:
        return coP + (x - x2)/(x - x1)*P

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
    return (coP + 1/k*P)*M*(coP + k*P) - d/k*P

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
    assuming Poincare rank of M at that point is p.

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
    assuming Poincare rank of M at that point is p.

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

def find_dual_basis_spanning_left_invariant_subspace(A, U):
    """Find matrix V, such that it's rows span a left invariant
    subspace of A and form a dual basis with columns of U (that
    is, V*U=I).
    """
    evlist = []
    for eigenval, eigenvects, evmult in A.eigenvectors_left():
        evlist.extend(eigenvects)
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
    combinedT = identity_matrix(M.base_ring(), M.nrows())
    poincare_map = singularities(M, x)
    reduction_points = [pt for pt,p in poincare_map.iteritems() if p >= 1]
    reduction_points.sort()
    singular_points = poincare_map.keys()
    singular_points.sort()
    while reduction_points:
        pointidx = rng.randint(0, len(reduction_points) - 1)
        point = reduction_points[pointidx]
        prank = poincare_map[point]
        while True:
            A0 = matrix_c0(M, x, point, prank)
            if A0.is_zero(): break
            A1 = matrix_c1(M, x, point, prank)
            U, V = alg1x(A0, A1, x)
            rng.shuffle(singular_points)
            for p2 in singular_points:
                if p2 == point: continue
                prank2 = poincare_map[p2]
                B0 = matrix_c0(M, x, p2, prank2)
                assert not B0.is_zero()
                v = find_dual_basis_spanning_left_invariant_subspace(B0, U)
                if v is not None:
                    point2 = p2
                    P = U*v
                    break
            else:
                point2 = any_integer(rng, M.base_ring(), poincare_map)
                P = U*V
            P = P.simplify_rational()
            M = balance_transform(M, P, point, point2, x)
            M = M.simplify_rational()
            combinedT = combinedT * balance(P, point, point2, x)
            if point2 not in poincare_map:
                poincare_map[point2] = 1
                singular_points.append(point2)
        poincare_map[point] = prank - 1
        if prank <= 1:
            del reduction_points[pointidx]
    combinedT = combinedT.simplify_rational()
    return M, combinedT

def is_singularity_apparent(M0):
    """Return True if a matrix residue M0 corresponds to an
    apparent singularity, False otherwise.

    The actual test checks if M0 is a diagonalizable matrix with
    integer eigenvalues.
    """
    for eigenval, eigenvects, evmult in M0.eigenvectors_right():
        if len(eigenvects) != evmult:
            return False
        if not eigenval.is_integer():
            return False
    return True

def _find_balance_helper(a0, b0, x, x1, x2, cond1, cond2):

    def find_eval_by_cond(a0_ev, b0_ev, cond):
        for a0_eval, a0_evecs, a0_evmult in a0_ev:
            for b0_eval, b0_evecs, b0_evmult in b0_ev:
                if not cond(a0_eval, b0_eval):
                    continue
                for a0_evec in a0_evecs:
                    for b0_evec in b0_evecs:
                        if dot_product(a0_evec, b0_evec) == 1:
                            return a0_eval, a0_evec, b0_eval, b0_evec 
        return None, None, None, None

    a0_evl, b0_evr = a0.eigenvectors_left(), b0.eigenvectors_right()
    a0_eval, a0_evec, b0_eval, b0_evec = find_eval_by_cond(a0_evl, b0_evr, cond1)

    if (a0_evec and b0_evec) is None:
        a0_evr, b0_evl = a0.eigenvectors_right(), b0.eigenvectors_left()
        a0_eval, a0_evec, b0_eval, b0_evec = find_eval_by_cond(a0_evr, b0_evl, cond2)

        if (a0_evec and b0_evec) is None:
            return None

    P = cross_product(a0_evec, b0_evec)
    if cond1(a0_eval, b0_eval):
        T0 = balance(P, x2, x1, x)
    elif cond2(a0_eval, b0_eval):
        T0 = balance(P, x1, x2, x)

    return T0

def find_balance(a0, b0, x, x1, x2):
    """Check if the Fuchsian point x1 can be balanced with the singular point x2.
    If so, return the transformatin T0, othervise return None.

    The algorithm is described in [1, Definition 5].
    """
    cond1 = lambda a0_eval, b0_eval: b0_eval < -0.5
    cond2 = lambda a0_eval, b0_eval: b0_eval >= 0.5

    T0 = _find_balance_helper(a0, b0, x, x1, x2, cond1, cond2)
    return T0


def find_mutual_balance(a0, b0, x, x1, x2):
    """Check if two Fuchsian points x1 and x2 can be mutually balanced.
    If so, return the transformatin T0, othervise return None.

    The algorithm is described in [1, Definition 6].
    """
    cond1 = lambda a0_eval, b0_eval: a0_eval >= 0.5 and b0_eval < 0.5
    cond2 = lambda a0_eval, b0_eval: a0_eval < 0.5 and b0_eval >= 0.5

    T0 = _find_balance_helper(a0, b0, x, x1, x2, cond1, cond2)
    return T0

def normalize(m, x):

    def print_eigenvalues(m, points):
        for xi in points:
            mi = matrix_c0(m, x, xi, 0)
            logger.info("  x = %s:" % xi)
            logger.info("      %s" % str(mi.eigenvalues()).replace("\n"," "))
            #logger.info("      %s" % str(mi.eigenvectors_left()).replace("\n"," "))
            #logger.info("      %s" % str(mi.eigenvectors_right()).replace("\n"," "))

    logger.info("The initial matrix:\n%s" % str(m))

    fuchs_points = [
            pt for pt,p in singularities(m, x).iteritems()
            if not is_singularity_apparent(matrix_c0(m, x, pt, p))
    ]
    T = identity_matrix(m.base_ring(), m.nrows())

    for x0 in fuchs_points:

        points = singularities(m, x).keys()
        apparent_points = [p for p in points if p not in fuchs_points]
        print points, fuchs_points, apparent_points
        for x2 in apparent_points:
            while True:
                logger.info("Balance a singular point x = %s with a Fuchsian point x = %s" % (x2, x0))
    
                a0 = matrix_c0(m, x, x0, 0)
                b0 = matrix_c0(m, x, x2, 0)
                logger.info("Eigenvalues:")
                points = singularities(m, x).keys()
                print_eigenvalues(m, points)
    
                T0 = find_balance(a0, b0, x, x0, x2)
                if T0 is None:
                    logger.info("No balance found!")
                    break
                else:
                    m = transform(m, x, T0)
                    T = T * T0
                    points = singularities(m, x).keys()
                    apparent_points = list(set(points) - set(fuchs_points))
                    if x2 not in apparent_points:
                        logger.info("A singular point x = %s succesfully balanced!\n" % x2)
                        break
                    else:
                        logger.info("Balance found!\n")

    points = singularities(m, x).keys()
    for x1, x2 in combinations(fuchs_points, 2):
        logger.info("Mutually balance Fuchsian points x = %s and x = %s" % (x1, x2))
        while True:
            a0 = matrix_c0(m, x, x1, 0)
            b0 = matrix_c0(m, x, x2, 0)
            logger.info("Eigenvalues:")
            points = singularities(m, x).keys()
            print_eigenvalues(m, points)

            T0 = find_mutual_balance(a0, b0, x, x1, x2)
            if T0 is None:
                logger.info("No balance found!")
                break
            else:
                m = transform(m, x, T0)
                T = T * T0
                points = singularities(m, x).keys()
                apparent_points = list(set(points) - set(fuchs_points))
                if (x1 not in points) or (x2 not in points):
                    logger.info("Fuchsian points x = %s and x = %s succesfully balanced!\n" % (x1, x2))
                    logger.info("Remaining singular points:")
                    print_eigenvalues(m, points)
                    break
                else:
                    logger.info("Balance found!\n")

    logger.info("\nThe transformation matrix:\n%s" % str(T))

    return m, T

def factor_epsilon(M, x, epsilon, seed=0):
    """Given a normalized Fuchsian system of differential equations:
        dF/dx = M(x,epsilon)*F,
    try to find a transformation that will factor out an epsilon
    from M. Return a transformed M (proportional to epsilon) and T.
    """
    n = M.nrows()
    rng = Random(seed)
    mu = SR.symbol()
    T_symbols = [SR.symbol() for i in xrange(n*n)]
    T = matrix(SR, n, n, T_symbols)
    eqs = []
    for point, prank in singularities(M, x).iteritems():
        assert prank == 0
        R = matrix_c0(M, x, point, 0)
        eq = (R/epsilon)*T-T*(R.subs({epsilon: mu})/mu)
        eqs.extend(eq.list())
    for solution in solve(eqs, T_symbols, solution_dict=True):
        # What follows is a kludge, but SageMath (as of 7.0)
        # can't handle symbolic temporary variables properly.
        # Here's an example:
        #     sage: t = SR.symbol(); t
        #     symbol11660
        #     sage: sol = solve(t*2=4, t, solution_dict=True)[0]; sol
        #     {symbol11660: 2}
        #     sage: t in sol
        #     False
        T = matrix([[solution[SR.symbols[str(e)]] for e in row] for row in T])
        if not T.is_invertible():
            continue
        # Right now T likely has a number of free variables in
        # it, we can set them to arbibtrary values, as long as
        # it'll make T invertible.
        while True:
            try:
                sT = T.subs([
                    e==rng.randint(-99, 99)
                    for e in T.variables() if e != epsilon
                ])
                sT = sT.simplify_rational()
                M = transform(M, x, sT)
            except (ZeroDivisionError, ValueError):
                continue
            break
        # We're leaking a bunch of temprary variables here,
        # which accumulate in SR.variables, but do I care?
        return M, sT
        # No, I don't.
    raise ValueError("can not factor epsilon")
