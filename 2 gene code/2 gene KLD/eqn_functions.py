# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:17:00 2018

@author: John
"""

# This file contains functions related to the model's SDE, including the functions governing the
# deterministic and noise dynamics. 

# These functions depend on the parameters found in eqn_parameters.py . 

import numpy as np
from eqn_parameters import kp, G, dm, dp, Nx, Ny, k, c_xx, c_yy, c_xy, c_yx, f_xx, f_yx, f_yy, f_xy

# ==============================================================

# Hill functions from deterministic part of protein SDE
def hill_x(x):
    
    hill_top = 0
    hill_bot = 0
    for i in range(0, Nx + 1):
        for j in range(0, Ny + 1):
            hill_top = hill_top + k[0, i, j]*c_xx[i, j]*c_xy[j]*(x[0]**i)*(x[1]**j)
            hill_bot = hill_bot + c_xx[i, j]*c_xy[j]*(x[0]**i)*(x[1]**j)
    
    hill_func_x = ((kp[0]*G[0])/dm[0])*(hill_top/hill_bot)
    return hill_func_x

def hill_y(x):
    
    hill_top = 0
    hill_bot = 0
    for i in range(0, Nx + 1):
        for j in range(0, Ny + 1):
            hill_top = hill_top + k[1, i, j]*c_yy[i, j]*c_yx[i]*(x[0]**i)*(x[1]**j)
            hill_bot = hill_bot + c_yy[i, j]*c_yx[i]*(x[0]**i)*(x[1]**j)
    
    hill_func_y = ((kp[1]*G[1])/dm[1])*(hill_top/hill_bot)
    return hill_func_y




# Decay terms in protein SDE.
def decay_x(x):
    decay_x = dp[0]*x[0]
    return decay_x

def decay_y(x):
    decay_y = dp[1]*x[1]
    return decay_y


# Deterministic SDE functions
def f1(x, params):
    result = hill_x(x) - decay_x(x)
    return result

def f2(x, params):
    result = hill_y(x) - decay_y(x)
    return result






# Binding-related noise term with f_xx in it
def b_xx(x):
    hill_top = 0
    hill_bot = 0
    
    for i in range(0, Nx):
        for j in range(0, Ny):
            hill_top = hill_top + 2*G[0]*f_xx[i, j]*c_xx[i, j]*c_xy[j]*(x[0]**i)*(x[1]**j)
    
    for i in range(0, Nx + 1):
        for j in range(0, Ny + 1):
            hill_bot = hill_bot + c_xx[i, j]*c_xy[j]*(x[0]**i)*(x[1]**j)
    
    result = (hill_top/hill_bot)*x[0]
    return result

# Binding-related noise term with f_yx in it
def b_yx(x):
    hill_top = 0
    hill_bot = 0
    
    for i in range(0, Nx):
        for j in range(0, Ny):
            hill_top = hill_top + 2*G[1]*f_yx[i, j]*c_yy[i, j]*c_yx[i]*(x[0]**i)*(x[1]**j)
    
    for i in range(0, Nx + 1):
        for j in range(0, Ny + 1):
            hill_bot = hill_bot + c_yy[i, j]*c_yx[i]*(x[0]**i)*(x[1]**j)
    
    result = (hill_top/hill_bot)*x[0]
    return result

# Binding-related noise term with f_yy in it
def b_yy(x):
    hill_top = 0
    hill_bot = 0
    
    for i in range(0, Nx):
        for j in range(0, Ny):
            hill_top = hill_top + 2*G[1]*f_yy[i, j]*c_yy[i, j]*c_yx[i]*(x[0]**i)*(x[1]**j)
    
    for i in range(0, Nx + 1):
        for j in range(0, Ny + 1):
            hill_bot = hill_bot + c_yy[i, j]*c_yx[i]*(x[0]**i)*(x[1]**j)
    
    result = (hill_top/hill_bot)*x[1]
    return result

# Binding-related noise term with f_xy in it
def b_xy(x):
    hill_top = 0
    hill_bot = 0
    
    for i in range(0, Nx):
        for j in range(0, Ny):
            hill_top = hill_top + 2*G[0]*f_xy[i, j]*c_xx[i, j]*c_xy[j]*(x[0]**i)*(x[1]**j)
    
    for i in range(0, Nx + 1):
        for j in range(0, Ny + 1):
            hill_bot = hill_bot + c_xx[i, j]*c_xy[j]*(x[0]**i)*(x[1]**j)
    
    result = (hill_top/hill_bot)*x[1]
    return result





def g1_gill(x, params):
    
    inside = f1(x, params) + decay_x(x) + b_xx(x) + b_yx(x)
    result = np.sqrt(inside)
    return result

def g2_gill(x, params):
    
    inside = f2(x, params) + decay_y(x) + b_yy(x) + b_xy(x)
    result = np.sqrt(inside)
    return result    
    
# ==============================================================
    



# =============================================
    


# Additive noise functions
    
def g1_add(x, params):
    # Additive noise
    return params[0]

def g2_add(x, params):
    # Additive noise
    return params[1]






# Multiplicative noise functions
    
def g1_mult(x, params):
    return params[2]*x[0] 

def g2_mult(x, params):
    return params[3]*x[1]
