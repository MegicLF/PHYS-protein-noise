# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:36:48 2018

@author: John
"""
# This file takes an SDE (deterministic fn, noise fn, derivative of noise fn) as input, along with a specified
# array of points, and evaluates the normalized solution to the 1D steady state Fokker-Planck equation at those
# points. 


import numpy as np


# =========================================================================

# IMPORTANT SUPPLEMENTARY FUNCTIONS




# This is the integral that must be calculated to solve the 1D steady state Fokker-Planck equation.
def whole_func(f, g, g_prime, x, params):
    numerator = g(x, params)*g_prime(x, params) - f(x, params)
    denominator = (g(x, params)**2)/2
    
    result = numerator/denominator
    return result


# Trapezoid rule integration. 
def int_trapezoid(f, g, g_prime, xmax, naive_step, int_min, params):
    
    num_pts = int(np.ceil((xmax - int_min)/naive_step))
    
    x = np.linspace(int_min, xmax, num_pts)
    dx = x[1] - x[0]
    
    middle = 0
    for i in range(1, num_pts-1):
        middle = middle + whole_func(f, g, g_prime, x[i], params)
    
    integral = (dx/2)*((whole_func(f, g, g_prime, x[0], params)) + 2*middle + whole_func(f, g, g_prime, x[num_pts - 1], params))
    
    return integral


# Normalizes Fokker-Planck solution.
def pss_norm(p, dx):
    num_pts = len(p)
    
    middle = 0
    for i in range(1, num_pts-1):
        middle = middle + p[i]
    
    integral = (dx/2)*(p[0] + 2*middle + p[num_pts - 1])
    
    norm = integral
    return p/norm

# =================================================================================
    


# MAIN FOKKER-PLANCK SOLVING FUNCTION

# f : deterministic function from protein SDE
# noise : noise function from protein SDE
# noise_deriv : derivative of noise function from protein SDE

# num_x_pts : number of points FP solution will be evaluated at
# x_min : smallest point FP solution evaluated at
# x_max : largest point FP solution evaluated at
# naive_step : step size used in trapezoidal integration (possibly adjusted to be smaller)
# int_min : Lower bound of V(x) integral. Useful to have it be > 0 if g(x) goes to zero for x -> 0. 
    
# params : additional parameters. Used here to specify sig_a for + noise sims, and sig_m for x noise sims.
    
def FP_solve(f, noise, noise_deriv, num_x_pts, x_min, x_max, naive_step, int_min, params):
    
    # Set up arrays.
    x = np.linspace(x_min, x_max, num_x_pts)       # Holds values FP solution will be evaluated at
    dx = x[1] - x[0]                               # Used for normalizing P_ss via trapezoid integration later. 
    
    V = np.zeros(num_x_pts)                        # Holds 'potential' used in FP solution.
    Pss = np.zeros(num_x_pts)                      # Holds P_ss. 
    
    
    #========================================================================
    
    
    # Calculate V for each desired x.
    for i in range(0, num_x_pts):
        V[i] = int_trapezoid(f, noise, noise_deriv, x[i], naive_step, int_min, params)
    
    
    # Calculate naive P_ss from V.
    for i in range(0, num_x_pts):
        Pss[i] = np.exp(-V[i])
    
    
    
    # Normalize P_ss.
    Pss_normed = pss_norm(Pss, dx)

    return Pss_normed