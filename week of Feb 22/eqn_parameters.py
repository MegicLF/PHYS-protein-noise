# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:11:05 2018
@author: John
"""

# This file contains all model parameters. Change them here.

import numpy as np



# Equation parameters.
kp = 2                     # Protein translation rate
G = 1                      # Total number of available gene sites 
dm = 1                     # mRNA degradation rate
dp = 1                     # protein degradation rate 
#N = 6                      # order of feedback (how many times protein can bind to gene) 
Nmin = 0
Nmax = 1

kmin = 20                  # Smallest mRNA transcription rate
kmax = 300                 # Largest mRNA transcription rate 
k = np.linspace(kmin, kmax, Nmax+1)       # Assume mRNA transcription rate increases linearly 
                                       # with number of proteins bound to gene


c_val = 0.1                # Ratio: forward gene-protein binding rate / backward gene-protein binding rate
c = np.ones(Nmax+1)           
for i in range(1,Nmax+1):
    c[i] = c_val**i        # Each binding event is independent, and all forward/backward ratios are equal,
                           # causing c_i = (c)^i

b_r_val = 0.1              # Backward gene-protein binding rate
b_r = np.zeros(Nmax+1)        
for i in range(1,Nmax+1): 
    b_r[i] = b_r_val       # For simplicity, assume all backward gene-protein binding rates are equal. 