import delirium.expression as expression

def alphabet(m, x):
    result = set()
    for ex in m.list():
        pts = expression.alphabet(ex, x)
        result.update(pts)
    return result

def expand_to_laurent_series(m, x):
    return m

def singular_points(m, x):
    result = set()
    for ex in m.list():
        pts = expression.singular_points(ex, x)
        result.update(pts)
    return result
