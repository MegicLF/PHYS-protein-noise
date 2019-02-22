# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 20:01:56 2019
@author: John
"""

import numpy as np

# ====================================================
    
# The SDE being simulated is
# x1-dot = f1(x1, x2) + g1(x1, x2) eta(t)
# x2-dot = f2(x1, x2) + g2(x1, x2) eta(t)
    
# x: 2 component vector representing current system state.
# f: 2 component vector containing deterministic parts of SDE
# g: 2 component vector containing stochastic parts of SDE
# step_size: SDE integration step size
# params: parameters used by f, g
    
def euler_maruyama_2D(x, f, g, step_size, params):
    D = len(x)
    x_new = -100*np.ones(D)
    
    for i in range(0, D):      
        while x_new[i] < 0:
            r = np.random.normal(0,1)
            x_new[i] = x[i] + f[i](x, params)*step_size + g[i](x, params)*np.sqrt(step_size)*r
   
    return x_new  
