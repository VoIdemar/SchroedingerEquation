import types
import numpy as np

def simpson_method(f, a, b, N=None, eps=None, dN=4):
    """
    Numerically computes the value of the definite integral of the function f
    over the interval (a, b) ([a, b]) using Simpson's rule
    """
    def mesh_simpson(f, a, b, N):
        h = (b - a)/float(N)
        x = np.linspace(a, b, N + 1)
        
        g = ((lambda k: f(x[k])) if isinstance(f, types.FunctionType) else 
             (lambda k: f[k]))   
        
        integral = 0
        for k in range(1, N, 2):
            integral = integral + g(k-1) + 4*g(k) + g(k+1)
        return h*integral/3.0
        
    def eps_mesh_simpson(f, a, b, eps, dN, N0=10):
        N = N0 + dN
        i1 = mesh_simpson(f, a, b, N0)
        i2 = mesh_simpson(f, a, b, N)
        while np.abs(i1-i2) > eps:
            N = N + dN
            i1 = i2
            i2 = mesh_simpson(f, a, b, N)
        return (i1 + i2)/2.0
    
    if eps is None and N is None:
        raise ValueError("Wrong function arguments - one of the (N, eps) arguments should be specified")
    elif eps is None:
        if N <= 0 or (N % 2) <> 0:
            raise ValueError("Wrong function arguments - N != 2*m, m > 0")
        return mesh_simpson(f, a, b, N)
    elif N is None:
        if eps <= 0 or eps >= 1:
            raise ValueError("Wrong function arguments - eps <= 0 or eps >= 1")
        if dN < 1:
            raise ValueError("Wrong function arguments - dN < 1")      
        return eps_mesh_simpson(f, a, b, eps, dN)
    else:
        raise ValueError("Ambiguity - N and eps specified. Exactly one required")