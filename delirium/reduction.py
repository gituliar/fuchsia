#!/usr/bin/env python
from sage.all import *
from random import Random

def transform(M, x, T):
    """Given a system of differential equations dF/dx = M*F,
    and a transformation of base functions F = T*F', compute
    and return M', such that dF'/dx = M'*F'.

    Note: M' = inverse(T)*(M*T - dT/dx)
    """
    return T.inverse()*(M*T - derivative(T, x))

def balance(P, x1, x2, x):
    assert P.is_square()
    assert (P*P - P).is_zero()
    assert x1 != x2
    coP = identity_matrix(P.nrows()) - P
    if x1 == oo:
        return coP - (x - x2)*P
    elif x2 == oo:
        return coP - 1/(x - x1)*P
    else:
        return coP + (x - x2)/(x - x1)*P

def singularities(M, x):
    """Find values of x around which rational matrix M has
    a singularity; return a dictionary with {val: k} entries,
    such that:
        lim(M, x->val) -> C/(x-val)**k

    Note 1: use matrix_limit(M*(x-val)**k, x=val) to find C.

    Note 2: singularity at x->infinity is not reported.

    Example:
    >>> x, y = var("x y")
    >>> M = matrix([[1/x, 0], [1/(x+y)**3, 1/(x+y)]])
    >>> s = singularities(M, x)
    >>> from pprint import pprint; pprint(s)
    {0: 1, -y: 3}
    >>> matrix_limit(M*(x+y)**3, x=-y)
    [0 0]
    [1 0]
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
                result[val] = max(result[val], k)
            else:
                result[val] = k
    return result

def matrix_limit(M, **kwargs):
    """Same as limit(), but works on matrices.

    Example:
    >>> x = var('x')
    >>> matrix_limit(matrix([[sin(x)/x, x]]), x=0)
    [1 0]
    """
    return matrix([[limit(e, **kwargs) for e in row] for row in M])

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
    if not (A*C - B).is_zero():
        raise ValueError, "matrix equation has no solutions"
    return C

def solve_left_fixed(A, B):
    """As of SageMath 6.10, 'Matrix.solve_left' method uses a broken
    check for solution correctness; this function corrects that check.
    """
    C = A.solve_left(B, check=False)
    if not (C*A - B).is_zero():
        raise ValueError, "matrix equation has no solutions"
    return C

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
    with M being singular around x=v like so: lim(M, x=v)->C/(x-v)**p,
    try to find a transformation T, which will reduce p by 1 (but
    possibly introduce another singularity at x=v2). Return the
    transformed M and T.
    """
    assert M.is_square()
    assert p > 1
    n = M.nrows()
    combinedT = identity_matrix(n)
    while True:
        A0 = matrix_limit(M*(x-v)**p, x=v)
        if A0.is_zero(): break
        A1 = matrix_limit((M - A0*(x-v)**(-p))*(x-v)**(p-1), x=v)
        U, V = alg1x(A0, A1, x)
        P = U*V
        T = balance(P, v, v2, x)
        M = transform(M, x, T)
        M = M.simplify_rational()
        combinedT = combinedT * T
    combinedT = combinedT.simplify_rational()
    return M, combinedT

def any_integer(rng, ring, excluded):
    r = 2
    while True:
        p = ring(rng.randint(-r, r))
        if p not in excluded:
            return p
        r *= 2

def find_dual_basis_in_invariant_space(A, U):
    """Find matrix V, such that it's rows belong to a left
    invariant space of A and form a dual basis with columns of U
    (that is, V*U=I).
    """
    for eigenval, eigenvects, evmult in A.eigenvectors_left():
        W = matrix(eigenvects)
        try:
            M = solve_left_fixed(W*U, identity_matrix(U.ncols()))
            return M*W
        except ValueError:
            pass
    return None

def fuchsianize(M, x, seed=0):
    """Given a system of differential equations of the form dF/dx=M*F,
    try to find a transformation T, which will reduce M to Fuchsian
    form. Return the transformed M and T.

    Note that such transformations are not unique; you can obtain
    different ones by supplying different seeds.
    """
    assert M.is_square()
    rng = Random(seed)
    combinedT = identity_matrix(M.base_ring(), M.nrows())
    exponent_map = singularities(M, x)
    reduction_points = [p for p,e in exponent_map.iteritems() if e >= 2]
    reduction_points.sort()
    while reduction_points:
        pointidx = rng.randint(0, len(reduction_points) - 1)
        point = reduction_points[pointidx]
        exp = exponent_map[point]
        while True:
            A0 = matrix_limit(M*(x-point)**exp, x=point)
            if A0.is_zero(): break
            A1 = matrix_limit(
                    (M-A0*(x-point)**(-exp))*(x-point)**(exp-1), x=point)
            U, V = alg1x(A0, A1, x)
            for p2, exp2 in exponent_map.iteritems():
                if p2 == point: continue
                B0 = matrix_limit(M*(x-p2)**exp2, x=p2)
                assert not B0.is_zero()
                v = find_dual_basis_in_invariant_space(B0, U)
                if v is not None:
                    point2 = p2
                    P = U*v
                    break
            else:
                point2 = any_integer(rng, M.base_ring(), exponent_map)
                P = U*V
            T = balance(P, point, point2, x)
            M = transform(M, x, T)
            M = M.simplify_rational()
            combinedT = combinedT * T
            if point2 not in exponent_map:
                exponent_map[point2] = 1
        exponent_map[point] = exp - 1
        if exp <= 2:
            del reduction_points[pointidx]
    combinedT = combinedT.simplify_rational()
    return M, combinedT
