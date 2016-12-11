import numpy as np

def numerov_method(a, b, N, E, q, psi0, psi1, backward=False):
    h = (b - a)/float(N)
    c = (h**2) / 12.0
    x = np.linspace(a, b, N + 1)
    
    psi = np.zeros(N + 1)
    if not backward:
        psi[0] = psi0
        psi[1] = psi1
    else:
        psi[N] = psi0
        psi[N-1] = psi1
        
    indexes = range(2, N + 1) if not backward else range(N - 2, -1, -1)
    qk = lambda k: c*q(E, x[k])
    
    if not backward:
        for k in indexes:
            psi[k] = (2.0*(1.0-5.0*qk(k-1))*psi[k-1] - (1.0+qk(k-2))*psi[k-2])/(1.0+qk(k))
    else:
        for k in indexes:
            psi[k] = (2.0*(1.0-5.0*qk(k+1))*psi[k+1] - (1.0+qk(k+2))*psi[k+2])/(1.0+qk(k))
    
    return psi