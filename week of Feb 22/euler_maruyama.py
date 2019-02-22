# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 20:01:56 2019
@author: John
"""

import numpy as np
x=0

# Purpose: Steps SDE forward in time.
# The SDE being simulated is
# x-dot = f(x) + g(x) eta(t)

# x: Number representing current system state.
# f: Deterministic part of SDE.
# g: Stochastic part of SDE.
# step_size: smaller means more accurate but longer computation times
# params: vector containing any parameters that need to be passed to f or g.

def euler_maruyama_1D(x, f, g, step_size, params): 
    r = np.random.normal(0,1)    # Generate a normally distributed random number with mean 0 and variance 1.
    
    x_new = x + f(x, params)*step_size + g(x, params)*np.sqrt(step_size)*r
    
    return max(x_new,0)    # max is there to prevent concentrations from going below zero



# ====================================================
    
# REQUIREMENTS:
# - Must return a 2 component vector x_new = (x1_new, x2_new)
# - No component of x_new must be less than or equal to zero.


    
# The SDE being simulated is
# x1-dot = f1(x1, x2) + g1(x1, x2) eta(t)
# x2-dot = f2(x1, x2) + g2(x1, x2) eta(t)
    
# x: 2 component vector representing current system state.
# f: 2 component vector containing deterministic parts of SDE
# g: 2 component vector containing stochastic parts of SDE
# step_size: same as above
# params: same as above
    
def euler_maruyama_2D(x, f, g, step_size, params):
    x1, x2 = x[0], x[1]
    f1, f2 = f[0], f[1]
    g1, g2 = g[0], g[1]
    r = np.random.normal(0,1)

    x1_new = x1 + f1(x1, x2, params)*step_size + g1(x1, x2, params)*np.sqrt(step_size)*r
    x2_new = x2 + f2(x1, x2, params)*step_size + g2(x1, x2, params)*np.sqrt(step_size)*r
    
    x_new = [max(x1_new, 0), max(x2_new, 0)]

    return x_new  

