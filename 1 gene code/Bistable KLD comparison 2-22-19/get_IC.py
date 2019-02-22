# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 12:26:44 2018

@author: John
"""

from eqn_functions import f

# Classical Runge-Kutta RK4 solver.
def RK4(x_old, f, step_size, params):
    j1 = f(x_old, params)
    j2 = f(x_old + (step_size/2)*j1, params)  
    j3 = f(x_old + (step_size/2)*j2, params)
    j4 = f(x_old + (step_size)*j3, params)
    
    x_new = x_old + (step_size/6)*(j1 + 2*j2 + 2*j3 + j4)
     
    return x_new


def get_IC(test_x, step_size, params):

    flag = 0
    while flag == 0:
        sol = RK4(test_x, f, step_size, params) # Runge-Kutta method
        
        if test_x != sol:
            test_x = sol
        else:
            flag = 1
    x0 = test_x
    
    return x0
    