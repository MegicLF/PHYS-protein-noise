# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 22:15:09 2018

@author: Linghui
"""

import numpy as np
import histogram


# =======================================================
# 1. PARAMETERS

# Reaction propensity parameters
k_low = 1.13           # low transcription rate
k_high = 1.67
k_p =  0.6          # translation

b_min = 0.1
b_max = 0.6
 
#need to change to smaller numbers
d_m = 1/60         # mRNA degradation 0.0017 0.280
d_p = 1/150            # protein degradation 0.00028~0.0017(0.084~0.51)


# Initial conditions for  g0, mRNA, protein
gt = 1  # number of types of gene involved
N=3
initial_g0 = 3
initial_g1 = 0
initial_g2 = 0
initial_g3 = 0
initial_g4 = 0
initial_g5 = 0
initial_g6 = 0
initial_m = 0
#initial_p = (k_low*initial_g0*k_p)/(d_m*d_p)  # p0=pss=(km*g0*kp)/(dm*dp)
G = 1

# Other parameters
vol = 4000    # Cell volume; used to convert between number and concentration


# Timestep parameters
num_timesteps = 10000
time_limit = 5000
step_size = (time_limit)/(num_timesteps)
data_size = int(0.10*step_size)

B = 1000  # ratio of b_forward and b_backward

TYPE = 0 #type of noise: additive-1, multiplicative-2, Gillespie-3

b_forward = np.linspace(b_min,b_max,6)
b_backward = [b_forward[i]/B for i in range(6)]
k_trans = np.linspace(k_low,k_high,6)
c = np.zeros(6)
c[0] = b_forward[0]/b_backward[0]
for i in range (1,6):
    c[i] = (b_forward[i]/b_backward[i])*c[i-1]
    

Num_sim = 1 # number of simulations
noise_size = 1
p_AddNoise = np.linspace(0.2,10,noise_size)
p_MultiNoise = np.linspace(0.01,0.50,noise_size)

histogram
