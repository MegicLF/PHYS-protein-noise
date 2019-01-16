# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:17:00 2018

@author: John
"""

# This file contains functions related to the model's SDE, including the functions governing the
# deterministic and noise dynamics. 

# These functions depend on the parameters found in eqn_parameters.py . 

import numpy as np
from eqn_parameters import kp, G, dm, dp, k, c, b_r



# ==============================================================
    

# Hill function from protein SDE related to protein translation.    
def hill(x, params):
    N = int(params[3])
    
    hill_top = 0
    hill_bot = 0
    for i in range(0, N + 1):
        hill_top = hill_top + k[i]*c[i]*(x**i)
        hill_bot = hill_bot + c[i]*(x**i)
    
    hill_func = ((kp*G)/dm)*(hill_top/hill_bot) 

    
    result = hill_func 
    return result

# Decay term in protein SDE.
def decay(x):
    decay = dp*x
    
    result = decay
    return result




# Hill function related to gene-protein binding from the noise part of the protein SDE.
def hill_bind(x, params):
    N = int(params[3])
    
    hill_top = 0
    hill_bot = 0
    for i in range(1,N+1):
        hill_top = hill_top + b_r[i]*c[i]*(x**i)
    
    for i in range(0,N+1):
        hill_bot = hill_bot + c[i]*(x**i)
        
    hill_func = 2*G*(hill_top/hill_bot)
    
    result = hill_func
    return result





# Overall deterministic part of the protein SDE.
def f(x, params):
    
    result = hill(x, params) - decay(x)
    return result

# Overall noise part of the protein SDE.
def noise(x, params):
    result = np.sqrt( hill(x, params) + decay(x) + hill_bind(x, params) )
    return result

# ====================================================================
    

# Derivative of the above function hill(x). Used in Fokker-Planck solution.
def hill_deriv_1(x, params):
    N = int(params[3])
    
    f1 = 0
    f2 = 0
    g1 = 0
    g2 = 0
    
    for i in range(0,N+1):
        f1 = f1 + k[i]*c[i]*(x**i)
        g1 = g1 + c[i]*(x**i)
        
    for i in range(1,N+1):
        f2 = f2 + i*k[i]*c[i]*(x**(i-1))
        g2 = g2 + i*c[i]*(x**(i-1))
        
    deriv = (f2*g1 - f1*g2)/(g1**2)
    
    
    result = ((kp*G)/dm)*deriv
    return result


# Derivative of the above function hill_bind(x). Used in Fokker-Planck solution.
def hill_deriv_2(x, params):
    N = int(params[3])
    
    f1 = 0
    f2 = 0
    g1 = 0
    g2 = 0
    
    for i in range(0,N+1):
        f1 = f1 + b_r[i]*c[i]*(x**i)
        g1 = g1 + c[i]*(x**i)
        
    for i in range(1,N+1):
        f2 = f2 + i*b_r[i]*c[i]*(x**(i-1))
        g2 = g2 + i*c[i]*(x**(i-1))
        
    deriv = (f2*g1 - f1*g2)/(g1**2)
    
    result = 2*G*deriv
    return result


# Overall derivative of noise(x). Used in Fokker-Planck solution.
def noise_deriv(x, params):
    numerator = hill_deriv_1(x, params) + dp + hill_deriv_2(x, params)
    denominator = 2*np.sqrt(noise(x, params))
    
    result = numerator/denominator
    return result




# =============================================
    


# Additive noise functions
    
def g_add(x, params):
    # Additive noise
    return params[0]


def g_prime_add(x, params):
    return 0






# Multiplicative noise functions
    
def g_mult(x, params):
    return params[1]*x + params[2]

def g_prime_mult(x, params):
    return params[1]