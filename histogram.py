# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 09:35:59 2018

@author: fengl
"""

import QSS
import numpy as np
import matplotlib.pyplot as plt

N = 10  #number of simulations

p_A = [[0] for a in range(N)]
p_M = [[0] for a in range(N)]
p_G = [[0] for a in range(N)]

for i in range(N):
    p_A[i] = QSS.p_add
    p_M[i] = QSS.p_multi
    p_G[i] = QSS.p_Gillespie

plt.hist(p_A,bins='auto',range=(QSS.initial_p-100,QSS.initial_p+100))
plt.show()
plt.hist(p_M,bins='auto',range=(QSS.initial_p-100,QSS.initial_p+100))
plt.show()
plt.hist(p_G,bins='auto',range=(QSS.initial_p-100,QSS.initial_p+100))
plt.show()