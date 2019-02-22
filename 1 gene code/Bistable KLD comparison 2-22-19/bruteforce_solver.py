# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:12:14 2018

@author: John
"""


import numpy as np
from scipy import stats
from solution_options import tsize, step_size, start_index, num_sims


# ====================================================================

# Step SDE forward in time.
def euler_maruyama(x,f,g,step_size, params): 
    r = np.random.normal(0,1)
    x_new = x + f(x, params)*step_size + g(x, params)*np.sqrt(step_size)*r
    return max(x_new,0)



# ===================================================================
    
# MAIN BRUTE-FORCE SOLVING FUNCTION
    
# Does a bunch of simulations, then uses some of the simulation data (the time points which are
# past the time step threshold specified by 'start_index') to generate a kernel density estimate (KDE).
# Then that KDE is evaluated at specified values and returned. 
    
# f : deterministic function from protein SDE
# noise : noise function from protein SDE

# num_x_pts : number of points KDE solution will be evaluated at
# x_min : smallest point KDE solution evaluated at. Pick -1 to use lowest simulation value.
# x_max : largest point KDE solution evaluated at. Pick -1 to use highest simulation value.
    
def bruteforce_solve(f, noise, protein_IC, num_xpts, x_min, x_max, params):
    
    # Set up simulation.
    protein = np.zeros((num_sims, tsize + 1))   # Will hold all protein SDE simulation data.
    
    for i in range(0, num_sims):
        protein[i,0] = protein_IC               # Initialize each simulation.
    
    
    
    # Run simulation.
    for j in range(0, num_sims):
        for i in range(0, tsize):
            protein[j,i+1] = euler_maruyama(protein[j,i],f,noise,step_size,params)
    
    
    flat_protein = protein[:,start_index:].flatten()   # Put all simulation data in 1D array
                                                       # (with time points below threshold 'start_index' removed) 
                                                       
    kernel = stats.gaussian_kde(flat_protein)          # Generate kernel density estimate using flattened array
    
    
    
    
    # If x_min isn't already specified, use minimum of simulation values.
    if x_min == -1:
        p_min = np.min(flat_protein)
    else:
        p_min = x_min
    
    # If x_max isn't already specified, use maximum of simulation values.
    if x_max == -1:
        p_max = np.max(flat_protein)
    else:
        p_max = x_max
    
    
    
        
    values = np.linspace(p_min, p_max, num_xpts)       # Generate array of values pdf will be evaluated at 
    KDE_pdf = kernel(values)                           # Evaluate KDE pdf at those values
    
    return KDE_pdf, p_min, p_max