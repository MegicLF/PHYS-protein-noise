# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:11:05 2018

@author: John
"""

# This file contains all model parameters. Change them here.

import numpy as np



# 2D equation parameters
Nx = 6
Ny = 6

kp = np.zeros(2)
G = np.zeros(2)
dm = np.zeros(2)
dp = np.zeros(2)

kmin = 1
kmax = 100
k = np.zeros((2, Nx + 1, Ny + 1))

c = 1/10
c_xx = np.zeros((Nx + 1, Ny + 1))
c_yy = np.zeros((Nx + 1, Ny + 1))
c_xy = np.zeros(Ny + 1)
c_yx = np.zeros(Nx + 1)

f_xx = np.zeros((Nx, Ny))
f_yy = np.zeros((Nx, Ny))
f_xy = np.zeros((Nx, Ny))
f_yx = np.zeros((Nx, Ny))


for i in range(0, 2):
    kp[i] = 1
    G[i] = 1
    dm[i] = 1
    dp[i] = 1
    
   # for j in range(0, Nx + 1):
    #    for l in range(0, Ny + 1):
            #k[i, j, l] = kmin + ((kmax - kmin)/(Nx + Ny))*(j + l)
            
#            if (j == Nx) or (l == Ny):
#                k[i, j, l] = kmax
#            else:
#                k[i, j, l] = 0

#for j in range(0, Nx + 1):
#    for l in range(0, Ny + 1):
#        if (j == Nx) and (l < Ny):
#            k[0, j, l] = kmax
#            k[1, j, l] = 0
#        elif (j == Nx) and (l == Ny):
#            k[0, j, l] = (kmin + kmax)/2
#            k[1, j, l] = (kmin + kmax)/2
#        elif (j < Nx) and (l == Ny):
#            k[0, j, l] = 0
#            k[1, j, l] = kmax
#        else:
#            k[0, j, l] = kmin
#            k[1, j, l] = kmin
            
            
for j in range(0, Nx + 1):
    for l in range(0, Ny + 1):
        if l > 0:
            k[0, j, l] = 0
        elif j == Nx:
            k[0, j, l] = kmax
        else: 
            k[0, j, l] = kmin
   
        if j > 0:
            k[1, j, l] = 0
        elif l == Ny:
            k[1, j, l] = kmax
        else:
            k[1, j, l] = kmin
            
            
for j in range(0, Nx + 1):
    for l in range(0, Ny + 1):
        c_xx[j, l] = c**(j + l)
        c_yy[j, l] = c**(j + l)        
        

for j in range(0, Nx + 1):
    c_yx[j] = 1

for l in range(0, Ny + 1):
    c_xy[l] = 1
    
    
for j in range(0, Nx):
    for l in range(0, Ny):
        f_xx[j, l] = 0.0001
        f_yy[j, l] = 0.0001
        f_xy[j, l] = 0.0001
        f_yx[j, l] = 0.0001
    
    