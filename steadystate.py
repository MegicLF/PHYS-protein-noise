# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 21:33:15 2018

@author: Linghui
"""

# find the steady state of p (under Gillespie noise)
import math
import parameter

k_p = parameter.k_p
G = parameter.G
d_m = parameter.d_m
d_p = parameter.d_p
gt = parameter.gt
k_min = parameter.k_low
k_max = parameter.k_high

k_trans = parameter.k_trans
c = parameter.c

step_size = parameter.step_size

lower_bound = (k_p/d_p)*(G/d_m)*k_min      # Lower bound guess for x_ss
upper_bound = (k_p/d_p)*(G/d_m)*k_max      # upper bound guess for x_ss

epsilon = 0.01

def my_func(p):
    part1 = 0.0
    part2 = 0.0
    for i in range(gt):
            part1 = part1 + k_trans[i]*c[i]*(p**i)
            part2 = part2 + c[i]*(p**i)        
    f = ((k_p*G)/d_m)*(part1/part2) - d_p*p
    return f


# Solve ODEs with fourth order (dt^4) accuracy
def RK4(p_old, f, step_size):
    j1 = f(p_old)
    
    j2 = f(p_old + (step_size/2)*j1)
    
    j3 = f(p_old + (step_size/2)*j2)
    
    j4 = f(p_old + (step_size)*j3)
    
    p_new = p_old + (step_size/6)*(j1 + 2*j2 + 2*j3 + j4)
    
    return p_new


# Function we use to get our x_ss. 
# Works like this: start at some random initial point test_x. 
    # Then keep timestepping that point until it doesn't change anymore.
    # When it stops changing, we've found our steady state.
def get_IC(test_p):

    p_old = test_p
    flag = 0
    
    while flag == 0:
        p_new = RK4(p_old, my_func, step_size) # Runge-Kutta method
        
        if (p_old != p_new):
            # If x_old doesn't equal x_new, keep going.
            p_old = p_new
        else:
            # If x_old = x_new, we've found the steady state!
            flag = 1
            
    p_ss = p_new
    
    return p_ss


print("Let's look for the steady state!\n")


my_ss1 = get_IC(lower_bound)
my_ss2 = get_IC(upper_bound)

N = round(math.log2((my_ss2-my_ss1)/epsilon))
my_ss3 = (my_ss2-my_ss1)/(2**N)


print("When I started at the lower bound p = "+str(lower_bound)+", I found:\np_ss = "+str(my_ss1))

print("When I started at the upper bound p = "+str(upper_bound)+", I found:\np_ss = "+str(my_ss2))

