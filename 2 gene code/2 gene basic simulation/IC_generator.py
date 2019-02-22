# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 18:42:26 2018

@author: John
"""
import numpy as np
from eqn_parameters import kmin, kmax



# Classical Runge-Kutta RK4 solver.
def RK4(x_old, f, step_size, param):
    j1 = np.zeros(2)
    j2 = np.zeros(2)
    j3 = np.zeros(2)
    j4 = np.zeros(2)
    
    for i in range(0, 2):
        j1[i] = f[i](x_old, param)
    
    for i in range(0, 2):
        j2[i] = f[i](x_old + (step_size/2)*j1, param)
    
    for i in range(0, 2):
        j3[i] = f[i](x_old + (step_size/2)*j2, param)
        
    for i in range(0, 2):
        j4[i] = f[i](x_old + (step_size)*j3, param)
        
  
    x_new = x_old + (step_size/6)*(j1 + 2*j2 + 2*j3 + j4)
  
    return x_new



def dist(x1, y1, x2, y2):
    t1 = (x1 - x2)**2
    t2 = (y1 - y2)**2
    
    dist = np.sqrt(t1 + t2)
    return dist

# =======================================================
    
def get_IC(f, test_x1, test_x2, step_size, stop_tol, param):

    test_x = [test_x1, test_x2]
    flag = 0
    while flag == 0:
        x_new = RK4(test_x, f, step_size, param) # Runge-Kutta method        
        
        if (test_x[0] == x_new[0]) and (test_x[1] == x_new[1]):
            flag = 1
        elif (dist(test_x[0], x_new[0], test_x[1], x_new[1]) < stop_tol):
            flag = 1
        else:
            test_x[0] = x_new[0]
            test_x[1] = x_new[1]
   
    x0 = test_x[0]
    y0 = test_x[1]
    
    return x0, y0
    



def generate_IC(tol, f, step_size, stop_tol, param):
    # tol used to address numerical error in finding attractors (same attractor may be found multiple times w/
    # slightly different coordinates)

    # step size used for finding attractor coordinates via RK4 integration of DEs
    
    c = param[4]
    
    
    # Look for 4 attractors.
    
    # x low, y high attractor
    #test_x_A = c/4
    #test_y_A = 2*c
    test_x_A = kmin
    test_y_A = kmin
    x_A, y_A = get_IC(f, test_x_A, test_y_A, step_size, stop_tol, param)
    print(x_A, y_A)
    
    # x high, y high attractor
    #test_x_B = 2*c
    #test_y_B = 2.1*c
    test_x_B = kmin
    test_y_B = kmax
    x_B, y_B = get_IC(f, test_x_B, test_y_B, step_size, stop_tol, param)
    print(x_B, y_B)
    
    # x low, y low attractor
    #test_x_C = c/4
    #test_y_C = c/5
    test_x_C = kmax
    test_y_C = kmin
    x_C, y_C = get_IC(f, test_x_C, test_y_C, step_size, stop_tol, param)
    print(x_C, y_C)
    
    # x high, y low attractor
    #test_x_D = 2*c
    #test_y_D = c/4
    test_x_D = kmax
    test_y_D = kmax
    x_D, y_D = get_IC(f, test_x_D, test_y_D, step_size, stop_tol, param)
    print(x_D, y_D)
    
    
    
    
    # See if any are basically the same.
    list_x = []
    list_y = []
    
    num_zeros = 1
    list_x.append(x_A)
    list_y.append(y_A)
    
    if dist(x_A, y_A, x_B, y_B) > tol:
        num_zeros = num_zeros + 1
        list_x.append(x_B)
        list_y.append(y_B)
        
        if (dist(x_A, y_A, x_C, y_C) > tol) and (dist(x_B, y_B, x_C, y_C) > tol):
            num_zeros = num_zeros + 1  
            list_x.append(x_C)
            list_y.append(y_C)
            
            if (dist(x_A, y_A, x_D, y_D) > tol) and (dist(x_B, y_B, x_D, y_D) > tol) and (dist(x_C, y_C, x_D, y_D) > tol):
                num_zeros = num_zeros + 1
                list_x.append(x_D)
                list_y.append(y_D)
        else:
            if (dist(x_A, y_A, x_D, y_D) > tol) and (dist(x_B, y_B, x_D, y_D) > tol):
                num_zeros = num_zeros + 1
                list_x.append(x_D)
                list_y.append(y_D)
            
    else: 
        if (dist(x_A, y_A, x_C, y_C) > tol):
                num_zeros = num_zeros + 1  
                list_x.append(x_C)
                list_y.append(y_C)  
                
                if  (dist(x_A, y_A, x_D, y_D) > tol) and (dist(x_C, y_C, x_D, y_D) > tol):
                    num_zeros = num_zeros + 1
                    list_x.append(x_D)
                    list_y.append(y_D)
        else:
            if dist(x_A, y_A, x_D, y_D) > tol:
                num_zeros = num_zeros + 1  
                list_x.append(x_D)
                list_y.append(y_D)
    
    
    protein_IC = []
    for i in range(0, num_zeros):
        protein_IC.append([list_x[i], list_y[i]])
    
    return num_zeros, protein_IC