def bisection_method(f, a, b, eps=0.001, max_iter=100):
    if a >= b:
        raise ValueError("a should be less than b")
    if f(a)*f(b) > 0:
        raise ValueError("Condition sgn(f(a)) != sgn(f(b)) not matched")
    n = 0
    L, R = a, b
    c = 0.0
    while n <= max_iter:
        n = n + 1
        c = (L + R)/2.0        
        fc = f(c)
        if fc == 0.0 or (R - L)/2.0 < eps:
            return c
        if fc*f(L) > 0:
            L = c
        else:
            R = c
    print "Maximum iteration reached"
    return c