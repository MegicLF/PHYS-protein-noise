# -*- coding: utf-8 -*-
"""
Created on Fri Oct 26 09:21:29 2018

@author: Linghui
"""

import pandas as pd
#import parameter
import math
import matplotlib.pyplot as plt

df1 = pd.read_csv('A=0.2M=0.01,6102.0 (0).csv')
df2 = pd.read_csv('A=10.0M=0.5,6102.0 (0).csv')

p1 = df1['protein(Gillespie noise)']
q1 = df1['protein(additive noise)']
p2 = df2['protein(Gillespie noise)']
q2 = df2['protein(additive noise)']
t = df1['time']

noise_size = 2
x = np.linspace(0.2,10,noise_size)

D1 = 0
D2 = 0
for i in range(len(t)):
    D1 = D1 + p1[i]*(math.log2(p1[i]/q1[i]))
    D2 = D2 + p2[i]*(math.log2(p2[i]/q2[i]))

D = [D1]+[D2]
    
plt.plot(x,D,color='blue')