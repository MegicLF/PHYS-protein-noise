# -*- coding: utf-8 -*-
"""
Created on Mon Sep 10 20:12:49 2018

@author: Linghui
"""

import numpy as np
import matplotlib.pyplot as plt
import random

# =======================================================
# 1. PARAMETERS

# Initial conditions for  g0, mRNA, protein
gt = 3  # number of types of gene involved
N=6
initial_g0 = 6
initial_g1 = 0
initial_g2 = 0
initial_g3 = 0
initial_g4 = 0
initial_g5 = 0
initial_g6 = 0
initial_m = 0
initial_p = 6102  # p0=pss=(km*g0*kp)/(dm*dp)
G = 1

# Reaction propensity parameters
k_low = 1.13           # low transcription rate
k_high = 1.67
k_p =  0.6          # translation
 
d_m = 1/60         # mRNA degradation 0.0017 0.280
d_p = 1/150            # protein degradation 0.00028~0.0017(0.084~0.51)

# Other parameters
vol = 4000    # Cell volume; used to convert between number and concentration


# Timestep parameters
num_timesteps = 10000
time_limit = 500
step_size = (time_limit)/(num_timesteps)


# =============================================================


# ===========================================================
# 2. Create and initialize t, g0, m, p arrays

g = [[0]*(num_timesteps + 1) for a in range(gt+1)]
m = np.zeros(num_timesteps + 1)
p = np.zeros(num_timesteps + 1)
#p_multi = np.zeros(num_timesteps + 1)
#p_Gillespie = np.zeros(num_timesteps + 1)
t = np.zeros(num_timesteps + 1)
k_trans = np.linspace(k_low,k_high,N+1)
b_forward = np.zeros(N)
b_backward = np.zeros(N)
c = np.zeros(N)

t[0] = 0

for b in range (gt):
    if b == 0:
        g[0][0] = initial_g0
    if b == 1:
        g[1][0] = initial_g1
    if b == 2:
        g[2][0] = initial_g2
    if b == 3:
        g[3][0] = initial_g3
    if b == 4:
        g[4][0] = initial_g4
    if b == 5:
        g[5][0] = initial_g5
    if b == 6:
        g[6][0] = initial_g6

m[0] = initial_m
p[0] = initial_p
#p_multi[0] = initial_p
#p_Gillespie[0] = initial_p
c[0] = 1



 
# ===================================================================

# =============================================================
# 3. METH-noiseODS-QSS+QSS
def QSS(p_noise):
    for i in range(0, num_timesteps):
        for j in range(gt+1):
           if j > 0:
               g[j][i+1] = updateGene(j,i)
           else:
               add = 0
               for k in range(1,gt):
                   add = add + g[k][i]
               g[0][i+1] = N-add
           if g[j][i+1] < 0:
               g[j][i+1] = 0
               
        m[i+1] = updateMRNA(N,i)
        p[i+1] = updateProtein(p[i],proteinDeterministic(N,c,i),proteinNoise(N,c,i),i,p_noise)
        t[i+1] = t[i] + step_size
        
        if m[i+1] < 0:
            m[i+1] = 0
        if p[i+1] < 0:
            p[i+1] = 0
        
    graph()
        
        
def updateGene(g_index,time):
    gene = (b_forward[g_index-1]/b_backward[g_index])*p[time]*g[g_index-1][time]
    return gene

def updateMRNA(N,time):
    Sum = 0.0
    for j in range(gt):
        Sum = Sum + k_trans[j]*g[j][time]
    return Sum/d_m
     
def updateProtein(p_old,deterministic,noise,time,p_noise):
    # use Euler-Maruyama method (SDE)
    # This simulates dy = f(y) dt + g(y) dW  for some function f. 
     
    eta = random.gauss(0,np.sqrt(step_size))
    # Get random number from Gauss dist with mean = 0 and var = step_size
    
    p_new = p_old + (step_size)*deterministic + noise*eta*p_noise
    #   p_(n+1) = p_n + deterministic_function(p_n) h + noise_function(p_n)*(random number)
    
    return p_new[0]

def proteinDeterministic(N,c,time):
    # part1 = k0 + k1*c1*p + k2*c2*p^2 +...+ kN*cN*p^N
    # part2 = 1 + c1*p + c2*p^2 +...+ cN*p^N
    # func = (kp*G/dm)*(part1/part2) - dp*p 
    part1 = 0.0
    part2 = 0.0
    for i in range(N):
        part1 = part1 + k_trans[i]*c[i]*(p**i)
        part2 = part2 + c[i]*(p**i)        
    func = ((k_p*G)/d_m)*(part1/part2) - d_p*p[time]
    return func

def proteinNoise(N,c,time):
    # part1 = k0 + k1*c1*p + k2*c2*p^2 +...+ kN*cN*p^N
    # part2 = 1 + c1*p + c2*p^2 +...+ cN*p^N
    # part3 = b10*c1*p + b21*c2*p^2 +...+ bN,N-1*cN*p^N
    # noise_function = sqrt((kp*G/dm)*(part1/part2) + dp*p + 2G*(part3/part4))
    part1 = 0.0
    part2 = 0.0
    part3 = 0.0
    for i in range(N):
        part1 = part1 + k_trans[i]*c[i]*(p**i)
        part2 = part2 + c[i]*(p**i)
        part3 = b_backward[i]*c[i]*(p**i)
    noise = ((k_p*G)/d_m)*(part1/part2) + d_p*p[time] + 2*G*(part3/part2)
    return np.sqrt(noise)
  
def graph():
    plt.plot(t,m,color='blue')
    plt.xlabel("Time")
    plt.ylabel("mRNA number")
    plt.title("mRNA number vs time")
    plt.show()
    
    plt.plot(t,p,color='red')
    plt.xlabel("Time")
    plt.ylabel("Protein concentration")
    plt.title("Protein concentration vs time")
    plt.show()

# =============================================================
b_forward = [0.2,0.4,0.3,0.2,0.1,0,0]
b_backward = [0,0.2,0.1,0.3,0.4,0.5,0.2]   
p_noise = 1
QSS(p_noise)