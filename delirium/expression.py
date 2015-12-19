from   sage.symbolic.ring import SR

def alphabet(ex, x):
    a, n = SR.wild(1), SR.wild(2)
    p_exp = a**n
    pts = singular_points(ex, x)
    result = set()
    for pt in pts:
        r = pt.match(p_exp)
        if r:
            result.add(r[a])
    return result

def singular_points(ex, x):
    #ex = ex.partial_fraction(x)
    pts = ex.denominator().operands()
    result = set()
    for pt in pts:
        if x not in pt.variables():
            continue
        result.add(pt)
    return result
