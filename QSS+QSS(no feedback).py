# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 12:09:50 2018

@author: Linghui
"""

import numpy as np
import matplotlib.pyplot as plt

# =======================================================
# 1. PARAMETERS

# Initial conditions for  g0, mRNA, protein
initial_g0 = 3
initial_m = 0
initial_p = 31
gt = 3

# Reaction propensity parameters
k_m0 = 1.13           # g0 transcription
k_p =  1.40          # translation
 
d_m = 0.51          # mRNA degradation 0.0017 0.280
d_p = 0.297            # protein degradation 0.00028~0.0017(0.084~0.51)

# Other parameters
vol = 4000    # Cell volume; used to convert between number and concentration


# Timestep parameters
num_timesteps = 100000
time_limit = 500

# for stat
T1=0.0
T2=0.0



# =============================================================


# ===========================================================
# 2. Create and initialize t, g0, m, p arrays

g01 = np.zeros(num_timesteps + 1)
m1 = np.zeros(num_timesteps + 1)
p1 = np.zeros(num_timesteps + 1)
t1 = np.zeros(num_timesteps + 1)


g02 = np.zeros(num_timesteps + 1)
m2 = np.zeros(num_timesteps + 1)
p2 = np.zeros(num_timesteps + 1)
t2 = np.zeros(num_timesteps + 1)


t1[0] = 0
g01[0] = initial_g0
m1[0] = initial_m
p1[0] = initial_p


t2[0] = 0
g02[0] = initial_g0
m2[0] = initial_m
p2[0] = initial_p

 
# ===================================================================


# =============================================================
# 3. METHODS-QSS+QSS
def QSS(T1,T2):
    for i in range(0, num_timesteps):
        # Calculate total propensity. 
        h01 = ( k_m0*g01[i]  + k_p*m1[i] + d_m*m1[i] + d_p*p1[i]  )
        #h02 = ( k_m0*g02[i] + k_m1*g12[i] + k_p*m2[i] + d_m*m2[i] + d_p*p2[i] 
                          #+ bf1*g02[i]*p2[i] + br1*g12[i] )
        
        # 1/h0 is the scale parameter.
        scale1 = 1/h01
        #scale2 = 1/h02
        
        # Randomly step forward in time.
        rand_tstep1 = 0
        #rand_tstep2 = 0
        while rand_tstep1==0:
            rand_tstep1 = np.random.exponential(scale1)
        t1[i+1] = rand_tstep1 + t1[i]
        #while rand_tstep2==0:
            #rand_tstep2 = np.random.exponential(scale2)
        t2[i+1] = rand_tstep1 + t2[i]
        
        if t1[i+1]<=time_limit:
            # g0 with constant 
            g01[i+1] = g01[i]
            
            m1[i+1] = k_m0*g01[i]/d_m
            
            p1[i+1] = protein_additive(g01[i],m1[i],p1[i],rand_tstep1)
         
        
            # Enforce positive numbers/concentrations.
        if t2[i+1]<=time_limit:
            g02[i+1] = g02[i]
            m2[i+1] = k_m0*g02[i]/d_m
            p2[i+1] = protein_Gillespie(g02[i],m2[i],p2[i],rand_tstep1)
            
            
        if g01[i+1] < 0:
            g01[i+1] = 0

        if m1[i+1] < 0:
            m1[i+1] = 0
        if p1[i+1] < 0:
            p1[i+1] = 0 
         
        
        if p2[i+1] < 0:
            p2[i+1] = 0
        if g02[i+1] < 0:
            g02[i+1] = 0
        if m2[i+1] < 0:
            m2[i+1] = 0

       
        if i>10:
            if (t1[i]-t1[i-10])/t1[i]<0.2 and T1==0:
                T1=t1[i]
                I=i
                
        if i>10:
            if (t2[i]-t2[i-10])/t2[i]<0.2 and T2==0:
                T2=t2[i]
                J=i
                
        
    graph()
    print("First timestep:")
    print(I)
    calculate1(I)
    calculate2(J)

    
def protein_additive(g0,m,p0,rand_tstep):
    # p SDE timestep
    p = (p0 + (k_p*m - d_p*p0)*rand_tstep + p_noise*np.random.normal(0,np.sqrt(rand_tstep)))
    return p


def protein_Gillespie(g0,m,p0,rand_tstep):
    # p SDE timestep
    p = (p0 + (k_p*m - d_p*p0)*rand_tstep + np.sqrt(k_p*m)*np.random.normal(0,np.sqrt(rand_tstep))-np.sqrt(d_p*p0)*np.random.normal(0,np.sqrt(rand_tstep)))
    return p

"""
def protein_multi(g0,m,p0,rand_tstep):
    # p SDE timestep
    p = (p0 + (k_p*m - d_p*p0)*rand_tstep + p_noise*p0*np.random.normal(0,np.sqrt(rand_tstep)))
    return p
"""

def calculate1(k):
    P1=np.zeros(num_timesteps-9+k+1)
    for j in range(k,num_timesteps-10):
        P1[j-k]=p1[j]
    mean=np.average(P1)  
    mean1.append(mean)
    SD=np.std(P1)
    SD1.append(SD)
    print("add_mean:")
    #print("multi_mean:")
    print(mean)
    print("add_SD:")
    #print("multi_SD:")
    print(SD)

def calculate2(k):
    P2=np.zeros(num_timesteps-9+k+1)
    for j in range(k,num_timesteps-10):
        P2[j-k]=p2[j]
    mean=np.average(P2) 
    mean2.append(mean)     
    SD=np.std(P2)
    SD2.append(SD)
    print("Gillespie_mean:")
    print(mean)
    print("Gillespie_SD:")
    print(SD)
            
# =============================================================
# 4. Plot results
def graph():
             
    plt.plot(t1,m1,color='blue')
    plt.xlabel("Time")
    plt.ylabel("mRNA number")
    plt.title("mRNA number vs time")
    plt.show()
    
    plt.plot(t1,p1,color='red')
    plt.plot(t2,p2,color='blue')
    plt.xlabel("Time")
    plt.ylabel("Protein concentration")
    plt.title("Protein concentration vs time")
    plt.show()

def graphStat():
    plt.plot(pNoise,mean1,color='red')
    plt.plot(pNoise,mean2,color='blue')
    plt.xlabel("p noise")
    plt.ylabel("mean")
    plt.title("p noise vs mean")
    plt.show()
    
    plt.plot(pNoise,SD1,color='red')
    plt.plot(pNoise,SD2,color='blue')
    plt.xlabel("p noise")
    plt.ylabel("SD")
    plt.title("p noise vs SD")
    plt.show()
#additive noise

mean1=[]
mean2=[]
SD1=[]
SD2=[]
pNoise=[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0]
p_noise=0.5
while p_noise<=10:
    QSS(T1,T2)
    print("p noise:")
    print(p_noise)
    p_noise=p_noise+0.5
graphStat()
"""
#multiplicative noise
p_noise=0.05
while p_noise<=0.20:
    QSS(T1,T2)
    print("p noise:")
    print(p_noise)
    p_noise=p_noise+0.01
"""

