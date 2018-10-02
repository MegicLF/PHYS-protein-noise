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
N=1
initial_g0 = 1
initial_g1 = 0
initial_g2 = 0
initial_g3 = 0
initial_g4 = 0
initial_g5 = 0
initial_g6 = 0
initial_m = 0
initial_p = (k_low*initial_g0*k_p)/(d_m*d_p)  # p0=pss=(km*g0*kp)/(dm*dp)
G = 1

# Other parameters
vol = 4000    # Cell volume; used to convert between number and concentration


# Timestep parameters
num_timesteps = 10000
time_limit = 500
step_size = (time_limit)/(num_timesteps)

B = 1000  # ratio of b_forward and b_backward

TYPE = 0 #type of noise: additive-1, multiplicative-2, Gillespie-3

# =============================================================


# ===========================================================
# 2. Create and initialize t, g0, m, p arrays

g_add = [[0]*(num_timesteps + 1) for a in range(gt+1)]
g_multi = [[0]*(num_timesteps + 1) for a in range(gt+1)]
g_Gillespie = [[0]*(num_timesteps + 1) for a in range(gt+1)]
m_add = np.zeros(num_timesteps + 1)
m_multi = np.zeros(num_timesteps + 1)
m_Gillespie = np.zeros(num_timesteps + 1)
p_add = np.zeros(num_timesteps + 1)
p_multi = np.zeros(num_timesteps + 1)
p_Gillespie = np.zeros(num_timesteps + 1)
t = np.zeros(num_timesteps + 1)
k_trans = np.linspace(k_low,k_high,N+1)
b_forward = np.linspace(b_min,b_max,6)
b_backward = [b_forward[i]/B for i in range(6)]
c = np.zeros(6)

t[0] = 0

for b in range (gt):
    if b == 0:
        g_add[0][0] = initial_g0
    if b == 1:
        g_add[1][0] = initial_g1
    if b == 2:
        g_add[2][0] = initial_g2
    if b == 3:
        g_add[3][0] = initial_g3
    if b == 4:
        g_add[4][0] = initial_g4
    if b == 5:
        g_add[5][0] = initial_g5
    if b == 6:
        g_add[6][0] = initial_g6

for b in range (gt):
    if b == 0:
        g_multi[0][0] = initial_g0
    if b == 1:
        g_multi[1][0] = initial_g1
    if b == 2:
        g_multi[2][0] = initial_g2
    if b == 3:
        g_multi[3][0] = initial_g3
    if b == 4:
        g_multi[4][0] = initial_g4
    if b == 5:
        g_multi[5][0] = initial_g5
    if b == 6:
        g_multi[6][0] = initial_g6
        
for b in range (gt):
    if b == 0:
        g_Gillespie[0][0] = initial_g0
    if b == 1:
        g_Gillespie[1][0] = initial_g1
    if b == 2:
        g_Gillespie[2][0] = initial_g2
    if b == 3:
        g_Gillespie[3][0] = initial_g3
    if b == 4:
        g_Gillespie[4][0] = initial_g4
    if b == 5:
        g_Gillespie[5][0] = initial_g5
    if b == 6:
        g_Gillespie[6][0] = initial_g6

m_add[0] = initial_m
m_multi[0] = initial_m
m_Gillespie[0] = initial_m

p_add[0] = initial_p
p_multi[0] = initial_p
p_Gillespie[0] = initial_p
c[0] = b_forward[0]/b_backward[0]
for i in range (1,6):
    c[i] = (b_forward[i]/b_backward[i])*c[i-1]


 
# ===================================================================

# =============================================================
# 3. METH-noiseODS-QSS+QSS
def QSS():
    for i in range(0, num_timesteps):
        updateAdd(i)
        updateMulti(i)
        updateGillespie(i)
        
        t[i+1] = t[i] + step_size
        
    graph(t,m_add,p_add)
    graph(t,m_multi,p_multi)
    graph(t,m_Gillespie,p_Gillespie)
        
def updateAdd(i):
        p_add[i+1] = updateProtein(p_add[i],proteinDeterministic(N,c,i,1),proteinAddNoise(N,c,i,p_AddNoise),i)
        
        m_add[i+1] = updateMRNA(N,i,1)
        
        if m_add[i+1] < 0:
            m_add[i+1] = 0
        if p_add[i+1] < 0:
            p_add[i+1] = 0
        
        for j in range(1, gt+1):
            g_add[j][i+1] = updateGene(j,i,p_add[i],g_add[j-1][i])
            if g_add[j][i+1] < 0:
                g_add[j][i+1] = 0
        S = 0 
        for k in range(gt):
            # S = ci*pi
            S = S + c[k]*p_add[k]
        g_add[0][i+1] = N/S
        if g_add[0][i+1] < 0: 
            g_add[0][i+1] = 0

def updateMulti(i):
        p_multi[i+1] = updateProtein(p_multi[i],proteinDeterministic(N,c,i,2),proteinMultiNoise(N,c,i,p_MultiNoise),i)
        
        m_multi[i+1] = updateMRNA(N,i,2)
        
        if m_multi[i+1] < 0:
            m_multi[i+1] = 0
        if p_multi[i+1] < 0:
            p_multi[i+1] = 0
        
        for j in range(1, gt+1):
            g_multi[j][i+1] = updateGene(j,i,p_multi[i],g_multi[j-1][i])
            add = 0     
            if g_multi[j][i+1] < 0:
                g_multi[j][i+1] = 0
        for k in range(gt):
            add = add + c[k]*p_multi[k]
        g_multi[0][i+1] = N/add
        if g_multi[0][i+1] < 0:
            g_multi[0][i+1] = 0

def updateGillespie(i):
        p_Gillespie[i+1] = updateProtein(p_Gillespie[i],proteinDeterministic(N,c,i,3),proteinGillNoise(N,c,i),i)
        
        m_Gillespie[i+1] = updateMRNA(N,i,3)
        
        if m_Gillespie[i+1] < 0:
            m_Gillespie[i+1] = 0
        if p_Gillespie[i+1] < 0:
            p_Gillespie[i+1] = 0
        
        for j in range(1, gt+1):
            g_Gillespie[j][i+1] = updateGene(j,i,p_Gillespie[i],g_Gillespie[j-1][i])
            add = 0     
            if g_Gillespie[j][i+1] < 0:
                g_Gillespie[j][i+1] = 0
        for k in range(gt):
            add = add + c[k]*p_Gillespie[k]
        g_Gillespie[0][i+1] = N/add
        if g_Gillespie[0][i+1] < 0:
            g_Gillespie[0][i+1] = 0


def updateGene(g_index,time,p,g):
    gene = (b_forward[g_index-1]/b_backward[g_index-1])*p*g
    return gene

def updateMRNA(N,time,TYPE):
    Sum = 0.0
    if TYPE == 1:
        for j in range(gt):
            Sum = Sum + k_trans[j]*g_add[j][time]
    if TYPE == 2:
        for j in range(gt):
            Sum = Sum + k_trans[j]*g_multi[j][time]
    if TYPE == 3:
        for j in range(gt):
            Sum = Sum + k_trans[j]*g_Gillespie[j][time]   
    return Sum/d_m
     
def updateProtein(p_old,deterministic,noise,time):
    # use Euler-Maruyama method (SDE)
    # This simulates dy = f(y) dt + g(y) dW  for some function f. 
     
    eta = random.gauss(0,np.sqrt(step_size))
    # Get random number from Gauss dist with mean = 0 and var = step_size
    
    p_new = p_old + (step_size)*deterministic + noise*eta
    #   p_(n+1) = p_n + deterministic_function(p_n) h + noise_function(p_n)*(random number)
    
    return p_new

def proteinDeterministic(N,c,time,TYPE):
    # part1 = k0 + k1*c1*p + k2*c2*p^2 +...+ kN*cN*p^N
    # part2 = 1 + c1*p + c2*p^2 +...+ cN*p^N
    # func = (kp*G/dm)*(part1/part2) - dp*p 
    part1 = 0.0
    part2 = 0.0
    if TYPE == 1:
        for i in range(gt):
            part1 = part1 + k_trans[i]*c[i]*(p_add[time]**i)
            part2 = part2 + c[i]*(p_add[time]**i)        
        func = ((k_p*G)/d_m)*(part1/part2) - d_p*p_add[time]
    if TYPE == 2:
        for i in range(gt):
            part1 = part1 + k_trans[i]*c[i]*(p_multi[time]**i)
            part2 = part2 + c[i]*(p_multi[time]**i)        
        func = ((k_p*G)/d_m)*(part1/part2) - d_p*p_multi[time]
    if TYPE == 3:
        for i in range(gt):
            part1 = part1 + k_trans[i]*c[i]*(p_Gillespie[time]**i)
            part2 = part2 + c[i]*(p_Gillespie[time]**i)        
        func = ((k_p*G)/d_m)*(part1/part2) - d_p*p_Gillespie[time]
    return func

def proteinAddNoise(N,c,time,p_AddNoise):
    return 1*p_AddNoise


def proteinMultiNoise(N,c,time,p_MultiNoise):
    return p_multi[time]*p_MultiNoise


def proteinGillNoise(N,c,time):
    # part1 = k0 + k1*c1*p + k2*c2*p^2 +...+ kN*cN*p^N
    # part2 = 1 + c1*p + c2*p^2 +...+ cN*p^N
    # part3 = b10*c1*p + b21*c2*p^2 +...+ bN,N-1*cN*p^N
    # noise_function = sqrt((kp*G/dm)*(part1/part2) + dp*p + 2G*(part3/part4))
    part1 = 0.0
    part2 = 0.0
    part3 = 0.0
    for i in range(gt):
        part1 = part1 + k_trans[i]*c[i]*(p_Gillespie[time]**i)
        part2 = part2 + c[i]*(p_Gillespie[time]**i)
        part3 = b_backward[i]*c[i]*(p_Gillespie[time]**i)
    noise = ((k_p*G)/d_m)*(part1/part2) + d_p*p_Gillespie[time] + 2*G*(part3/part2)
    return np.sqrt(noise)

  
def graph(t,m,p):
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

#def writeFile():
    
# =============================================================
p_AddNoise = 2
p_MultiNoise = 0.1
QSS()
