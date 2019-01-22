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
def euler_maruyama_2D(x,f,g,step_size, params): 
    
    
    # fill in your function here
    
    return x_new



# ===================================================================
    
# MAIN BRUTE-FORCE SOLVING FUNCTION
    
# Does a bunch of simulations, then uses some of the simulation data (the time points which are
# past the time step threshold specified by 'start_index') to generate a kernel density estimate (KDE).
# Then that KDE is evaluated at specified values and returned. 
    
# f = (f1, f2) : deterministic functions from protein SDE
# noise = (g1, g2) : noise functions from protein SDE

# num_x_pts : number of points KDE solution will be evaluated at (IN EACH DIMENSION)
# x_min = (x1_min, x2_min) : smallest point KDE solution evaluated at. Pick -1 to use lowest simulation value.
# x_max = (x1_max, x2_max) : largest point KDE solution evaluated at. Pick -1 to use highest simulation value.
    
def bruteforce_solve_2D(f, noise, protein_IC, num_xpts, x_min, x_max, params):
    
    # Set up simulation.
    protein = np.zeros((num_sims, tsize + 1, 2))   # Will hold all protein SDE simulation data.
    
    for i in range(0, num_sims):
        protein[i,0,0] = protein_IC[0]               # Initialize each simulation.
        protein[i,0,1] = protein_IC[1]
    
    
    
    # Run simulation.
    for j in range(0, num_sims):
        for i in range(0, tsize):
            protein[j,i+1,:] = euler_maruyama_2D(protein[j,i,:],f,noise,step_size,params)
    
    
    
    # FLAT PROTEIN: 2 components, flat_protein[0] holding all x simulation data, flat_protein[1] all y data
    flat_protein =                                    
                     # (with time points below threshold 'start_index' removed) 
                                                       
 
    
    
    p_min = np.zeros(2)
    p_max = np.zeros(2)
    # If x_min isn't already specified, use minimum of simulation values.
    for i in range(0,2):
        if x_min[i] == -1:
            p_min[i] = np.min(flat_protein[i])
        else:
            p_min[i] = x_min[i]
        
        # If x_max isn't already specified, use maximum of simulation values.
        if x_max[i] == -1:
            p_max[i] = np.max(flat_protein[i])
        else:
            p_max[i] = x_max[i]
    
    
    
    # Get KDE.
    X, Y = np.mgrid[p_min[0]:p_max[0]:np.complex(0, num_xpts), p_min[1]:p_max[1]:np.complex(0, num_xpts)]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([       ,            ])             #### PUT FLATTENED PROTEIN ARRAYS HERE
    kernel = stats.gaussian_kde(values)
    KDE_pdf = np.reshape(kernel(positions).T, X.shape)


    return KDE_pdf, p_min, p_max