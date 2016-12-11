"""
Numeric computation of the derivatives of the function f defined as a list of values
in the nodes of the mesh with step h in the node with index m
"""
def df_dx(f, m, h):
    return (f[m-2] - f[m+2] + 8.0*(f[m+1] - f[m-1]))/(12.0*h)

def df_dx_right(f, m, h):
    """ Computes right df_dx derivative """
    return (-11.0*f[m] + 18.0*f[m+1] - 9.0*f[m+2] + 2.0*f[m+3])/(6.0*h)

def df_dx_left(f, m, h):
    """ Computes left df_dx derivative """
    return (-2.0*f[m-3] + 9.0*f[m-2] - 18.0*f[m-1] + 11.0*f[m])/(6.0*h)

def d2f_dx2_center(f, m, h):
    return (f[m-1] - 2*f[m] + f[m+1])/(h**2)

def d2f_dx2_right(f, m, h):
    """ Computes right d2f_dx2 derivative """
    return (2*f[m] - 5*f[m+1] + 4*f[m+2] - f[m+3])/(h**2)

def d2f_dx2_left(f, m, h):
    """ Computes left d2f_dx2 derivative """
    return (-f[m-3] + 4*f[m-2] - 5*f[m-1] + 2*f[m])/(h**2)