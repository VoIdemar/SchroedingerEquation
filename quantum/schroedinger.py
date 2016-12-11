import numpy as np
import scipy.constants as const

from numerical.diffequ import numerov_method 
import numerical.differentiation as diff
from numerical.equ import bisection_method
from numerical.integration import simpson_method

def normalize_wave_func(wave_func, a, b, N):
    """
    Performs in-place quantum normalization of the wave function
    """
    c = simpson_method(wave_func**2, a, b, N)
    wave_func /= np.sqrt(c) 

def differentiate_wave_func(wave_func, a, b, N):
    h = (b - a)/float(N)
    return np.array([diff.df_dx_right(wave_func, 0, h), diff.df_dx_right(wave_func, 1, h)] +
                    [diff.df_dx(wave_func, k, h) for k in range(2, N - 1)] +
                    [diff.df_dx_left(wave_func, N - 1, h), diff.df_dx_left(wave_func, N, h)])

def find_radius_vector_expected_value(wave_func, a, b, N):
    x = np.linspace(a, b, N + 1)
    return simpson_method(wave_func*x*wave_func, a, b, N)

def find_delta_radius_vector(wave_func, a, b, N):
    x = np.linspace(a, b, N + 1)
    return x*wave_func - find_radius_vector_expected_value(wave_func, a, b, N)*wave_func

def find_quantum_expected_momentum(wave_func, a, b, N):
    """
    Calculates expectation value of the quantum momentum operator for the specified wave function   
    """ 
    dw_dx = differentiate_wave_func(wave_func, a, b, N)
    return -const.h*simpson_method(wave_func*dw_dx, a, b, N)

def find_quantum_expected_momentum_complex(wave_func, a, b, N):
    return complex(0.0, find_quantum_expected_momentum(wave_func, a, b, N))

def find_quantum_delta_momentum(wave_func, a, b, N):
    """
    Calculates quantum momentum delta operator for the specified wave function   
    """ 
    expected_momentum = find_quantum_expected_momentum(wave_func, a, b, N)
    momentum = -const.h*differentiate_wave_func(wave_func, a, b, N)
    return momentum - expected_momentum*wave_func

def find_quantum_delta_momentum_complex(wave_func, a, b, N):
    return complex(0.0, find_quantum_delta_momentum(wave_func, a, b, N))

def find_delta_expectation_values(wave_func, a, b, N):
    """
    Calculates expectation values of the radius vector and quantum momentum operators 
    for the specified wave function  
    """
    delta_rv = find_delta_radius_vector
    delta_momentum = find_quantum_delta_momentum
    quantum_dx_integral = delta_rv(delta_rv(wave_func, a, b, N), a, b, N)
    quantum_dx = simpson_method(wave_func*quantum_dx_integral, a, b, N)
    quantum_dmomentum_integral = delta_momentum(delta_momentum(wave_func, a, b, N), a, b, N)
    quantum_dmomentum = simpson_method(-wave_func*quantum_dmomentum_integral, a, b, N)
    return (quantum_dx, quantum_dmomentum)
    
def get_quantum_number(wave_func):
    """
    Gets quantum number of the specified wave function and the corresponding energy counting
    the number of sign changes
    """
    size = wave_func.shape[0]
    if size >= 2:
        qnum = 0
        # Boundary points are excluded to prevent wrong calculation of the quantum number
        for k in range(1, size - 2):
            s1, s2 = wave_func[k:k + 2]
            if s1*s2 < 0.0:
                qnum += 1
        return qnum
    else:
        return 0

def solve_schroedinger_eq(a, b, E, dE, U, d1, d2, N, match_node_idx, max_quantum_number):
    """
    Solves Schroedinger equation numerically using Numerov and bisection methods
    """
    tolerance = 0.0001
    h = (b - a)/float(N)
    q = lambda E, x: 2.0*(E - U(x))
    
    def is_correct_wave_func(wave_func, eps):
        return np.abs(wave_func[0]) <= eps and np.abs(wave_func[-1]) <= eps
    
    def solve_for_energy(E):
        """
        Finds the forward and the backward solutions of the Shroedinger equation and the difference
        between them in the match node
        """
        # Solving Shroedinger equation integrating forward
        wave_forward = numerov_method(a, b, N, E, q, 0.0, d1)
        # Solving Shroedinger equation integrating backward
        wave_backward = numerov_method(a, b, N, E, q, 0.0, d2, backward=True)
        
        # Common normalization
        factor = np.max(np.abs(wave_forward))
        wave_forward /= factor
        
        # Mathematical normalization        
        factor = wave_forward[match_node_idx]/wave_backward[match_node_idx]
        wave_backward *= factor
        
        return (wave_forward, wave_backward,
                diff.df_dx(wave_forward, match_node_idx, h) - diff.df_dx(wave_backward, match_node_idx, h))
    
    def find_energies_and_wave_functions(E0, dE, max_quantum_number, max_iter):
        f = lambda E: solve_for_energy(E)[2] 
        
        energies = []
        wave_functions = []
        qnumbers = []
        
        i = 0
        _, _, new_diff = solve_for_energy(E0)
        last_diff = 0
        Ei = E0
        while i <= max_iter:
            i = i + 1
            Ei = Ei + dE
            last_diff = new_diff
            _, _, new_diff = solve_for_energy(Ei)
            if last_diff*new_diff < 0:
                # Energy is specified by bisection method
                energy = bisection_method(f, Ei-dE, Ei, tolerance)
                # Corresponding wave function is calculated using Numerov method
                wave_func = numerov_method(a, b, N, energy, q, 0, d1)
                # Quantum normalization
                normalize_wave_func(wave_func, a, b, N)
                
                if not is_correct_wave_func(wave_func, tolerance):
                    continue
                
                # Calculating quantum number
                qnum = get_quantum_number(wave_func)
                if qnum > max_quantum_number:
                    break
                else:
                    print "New quantum level found:", energy
                    energies.append(energy)
                    wave_functions.append(wave_func)
                    qnumbers.append(qnum)
        _print_energy_levels(energies, qnumbers)
        return zip(qnumbers, energies, wave_functions)
                
    def _print_energy_levels(energies, quantum_numbers):
        for i in range(0, len(energies)):
            print str.format("Energy level #{0} = {1}", quantum_numbers[i], energies[i])   
         
    max_iter = int(100.0/dE)
    return find_energies_and_wave_functions(E, dE, max_quantum_number, max_iter)  