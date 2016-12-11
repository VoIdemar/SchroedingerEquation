import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt

from quantum.schroedinger import solve_schroedinger_eq, find_delta_expectation_values

def potential(x, v0, L, W):
    return v0*(W if abs(x) >= L else 
               -(x/float(L))**2)
    
def get_wave_func_by_qnum(qnum, shroedinger_res):
    for (qn, _, wave_func) in shroedinger_res:
        if qn == qnum:
            return wave_func
    return None

def plot_wave_functions(a, b, N, shroedinger_res, densities):
    
    def draw_subplot(subplot_id, title, x_axis_data, shroedinger_data):
        plt.subplot(subplot_id)
        plt.title(title)
        plt.axhline(linewidth=1, color='k', linestyle='--')
        plt.axvline(linewidth=1, color='k', linestyle='--')
        
        for (qnum, energy, wave_func) in shroedinger_data:
            plot_label = str.format("Energy level#{0} = {1}", qnum, energy)            
            plt.plot(x_axis_data, wave_func, label=plot_label)   
    
    x = np.linspace(a, b, N + 1)
    draw_subplot(211, "Wave functions in Schroedinger equation", x, shroedinger_res)
    draw_subplot(212, "Probability densities", x, densities)
    plt.subplots_adjust(left=0.04, right=0.8, top=0.96, bottom=0.04)
    plt.legend(bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, prop={'size':14})
    plt.show()

if __name__ == "__main__":
    N = 160
    match_node_idx = N/2 + 15
    L = 2.0
    d1 = 10**(-9)
    d2 = d1
    W = 3.0
    U0 = 1.0
    E0 = -1.0
    dE = 0.001
    max_quantum_number = 3
    #U = lambda x: U0 if abs(x) < L else W
    U = lambda x: potential(x, U0, L, W)
    results = solve_schroedinger_eq(-L, L, E0, dE, U, d1, d2, N, match_node_idx, max_quantum_number)
    densities = ((qnum, energy, wave_func**2) for (qnum, energy, wave_func) in results)
    wf_levels = {0 : get_wave_func_by_qnum(0, results), 
                 3 : get_wave_func_by_qnum(3, results)}
    for (qnum, wf_level) in wf_levels.items():
        if not (wf_level is None):
            x, p = find_delta_expectation_values(wf_level, -L, L, N)
            print "Expectation values level#", qnum, ":", (x, p)
            if x*p >= (const.h**2)/4.0:
                print "Uncertainty principle is satisfied for level", qnum
    plot_wave_functions(-L, L, N, results, densities)