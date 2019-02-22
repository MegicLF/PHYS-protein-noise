# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:36:09 2018
@author: John
"""

# Timer
import time as ti
start_time = ti.time()

import scipy.stats as stats
import numpy as np
import math

#import matplotlib as mpl           # uncomment these to run on ACCRE
#mpl.use('Agg')
import matplotlib.pyplot as plt

from solution_options import naive_step, int_min, num_KDE_pts, num_FP_pts, tsize, total_time, step_size
from eqn_functions import f, noise, noise_deriv, g_add, g_prime_add, g_mult, g_prime_mult
from get_IC import get_IC, RK4 
from eqn_parameters import kp, G, kmin, kmax, dp, dm, Nmin, Nmax


from FP_solver import FP_solve
from bruteforce_solver import bruteforce_solve,bruteforce_solve2

#from error_propagation import get_error

#import sys

# =====================================================================



N = np.arange(Nmin, Nmax + 1)
num_orders = len(N)


KDE_Gill = np.zeros((num_orders, num_KDE_pts))
KDE_Add = np.zeros((num_orders, num_KDE_pts))
KDE_Mult = np.zeros((num_orders, num_KDE_pts))
FP_Gill = np.zeros((num_orders, num_FP_pts))
FP_Add = np.zeros((num_orders, num_FP_pts))
FP_Mult = np.zeros((num_orders, num_FP_pts))


epsilon = 0
sig_num = 5
sig_a = np.linspace(10, 10, sig_num)
sig_m = np.linspace(0.01, 1, sig_num)


FP_KLD_Add = np.zeros(sig_num)
FP_KLD_Mult = np.zeros(sig_num)

KDE_KLD_Add = np.zeros(sig_num)
KDE_KLD_Mult = np.zeros(sig_num)



FP_KLD_Add_op = np.zeros(num_orders)
FP_KLD_Mult_op = np.zeros(num_orders)

KDE_KLD_Add_op = np.zeros(num_orders)
KDE_KLD_Mult_op = np.zeros(num_orders)


# IC lower/upper bounds and RK4 step size. 
pss_lowerbound = (kp*kmin*G)/(dp*dm)            # lower bound for p_ss
pss_upperbound = (kp*kmax*G)/(dp*dm)            # upper bound for p_ss
get_IC_stepsize = 0.01                          # step size used when finding IC  


for j in range(0, num_orders):
    for i in range(0, sig_num):
        # 1: sig_a, 2: sig_m, 3: epsilon, 4: order of feedback N
        params = [sig_a[i], sig_m[i], epsilon, N[j]]
    
    
        # -------------------------------------------------------------------------------------------------           
        
        pss_low = get_IC(pss_lowerbound, get_IC_stepsize, params)   # find closest steady state to lower bound
        pss_high = get_IC(pss_upperbound, get_IC_stepsize, params)  # find closest steady state to upper bound
        
        
        
        
        discrimination = 0.1                            # If |pss_low - pss_high| < this number, we say only 1 steady state.
                                                        # If |pss_low - pss_high| > this number, we say 2 steady states.
        
        difference = np.abs(pss_low - pss_high)
        if difference < discrimination:
            # Only 1 steady state. Difference between pss_low and pss_high probably due to numerical error.
            num_steadystates = 1   
            protein_IC = (pss_low + pss_high)/2
            protein_IC2 = protein_IC/4
        else:
            # Probably 2 steady states.
            num_steadystates = 2    
            protein_IC = [pss_low, pss_high]
            protein_IC2 = protein_IC/4
            
        print(protein_IC)
            
        # -------------------------------------------------------------------------------------------------    
        """
        KDE_Gill[j,:], p_min, p_max = bruteforce_solve(f, noise, protein_IC, num_KDE_pts, -1, -1, params)
        FP_Gill[j,:] = FP_solve(f, noise, noise_deriv, num_FP_pts, p_min, p_max, naive_step, int_min, params)

        KDE_Add[j,:], trash1, trash2 = bruteforce_solve(f, g_add, protein_IC, num_KDE_pts, p_min, p_max, params)
        FP_Add[j,:] = FP_solve(f, g_add, g_prime_add, num_FP_pts, p_min, p_max, naive_step, int_min, params)
    
        KDE_Mult[j,:], trash1, trash2 = bruteforce_solve(f, g_mult, protein_IC, num_KDE_pts, p_min, p_max, params)
        FP_Mult[j,:] = FP_solve(f, g_mult, g_prime_mult, num_FP_pts, p_min, p_max, naive_step, int_min, params)
        """
        
    # ============================================================================================
""" 
        for k in range(0, sig_num):
            FP_KLD_Add[k] = stats.entropy(FP_Gill[k,:], FP_Add[k,:])   
            #KDE_KLD_Add[k] = stats.entropy(KDE_Gill[k,:], KDE_Add[k,:])   
            
        for k in range(0, sig_num):
            FP_KLD_Mult[k] = stats.entropy(FP_Gill[k,:], FP_Mult[k,:])
            #KDE_KLD_Mult[k] = stats.entropy(KDE_Gill[k,:], KDE_Mult[k,:])
        
    FP_KLD_Add_op[j] = FP_KLD_Add.min() 
    #KDE_KLD_Add_op[j] = KDE_KLD_Add.min()
    FP_KLD_Mult_op[j] = FP_KLD_Mult.min()
    #KDE_KLD_Mult_op[j] = KDE_KLD_Mult.min()
    
    
vals_KDE = np.linspace(p_min, p_max, num_KDE_pts) # set up distribution
vals_FP = np.linspace(p_min, p_max, num_FP_pts)
"""

# =======================================================================

#plt.plot(vals_KDE, KDE_Gill, color='red')
#
#
#plt.plot(vals_FP, FP_Gill, color='blue', linewidth=3, linestyle='dashed')
#
#plt.title("KDE vs FP Gillespie comparison")
#plt.xlabel("x",fontsize=14)
#plt.ylabel("$P_{ss}$",fontsize=14)
#plt.savefig("N_KDE_vs_FP_Gill.png")
#plt.show()
#
#
#
#plt.plot(vals_FP, FP_Gill, color='black', linewidth=3, linestyle='dashed')
#for i in range(0, num_orders):
#    plt.plot(vals_FP, FP_Add[i,:], alpha = 0.5)
#plt.title("Add. noise comparison, FP")
#plt.xlabel("$x$",fontsize=14)
#plt.ylabel("$P_{ss}$",fontsize=14)
#plt.savefig("N_FP_Pss_Add.png")
#plt.show()
#
#plt.plot(vals_FP, FP_Gill, color='black', linewidth=3, linestyle='dashed')
#for i in range(0, num_orders):
#    plt.plot(vals_FP, FP_Mult[i,:], alpha = 0.5)
#plt.title("Mult. noise comparison, FP")
#plt.xlabel("$x$",fontsize=14)
#plt.ylabel("$P_{ss}$",fontsize=14)
#plt.savefig("N_FP_Pss_Mult.png")
#plt.show()
#
#
#
#
#
#plt.plot(vals_KDE, KDE_Gill, color='black', linewidth=3, linestyle='dashed')
#for i in range(0, num_orders):
#    plt.plot(vals_KDE, KDE_Add[i,:])
#plt.title("Add. noise comparison, KDE")
#plt.xlabel("$x$",fontsize=14)
#plt.ylabel("$P_{ss}$",fontsize=14)
#plt.savefig("N_KDE_Pss_Add.png")
#plt.show()
#
#plt.plot(vals_KDE, KDE_Gill, color='black', linewidth=3, linestyle='dashed')
#for i in range(0, num_orders):
#    plt.plot(vals_KDE, KDE_Mult[i,:])
#plt.title("Mult. noise comparison, KDE")
#plt.xlabel("$x$",fontsize=14)
#plt.ylabel("$P_{ss}$",fontsize=14)
#plt.savefig("N_KDE_Pss_Mult.png")
#plt.show()



# ====================================================================

"""
FP_KLD_Add = np.zeros(num_orders)
FP_KLD_Mult = np.zeros(num_orders)

KDE_KLD_Add = np.zeros(num_orders)
KDE_KLD_Mult = np.zeros(num_orders)

KLD_error_Add = np.zeros(num_orders)
KLD_error_Mult = np.zeros(num_orders)


for i in range(0, num_orders):
    FP_KLD_Add[i] = stats.entropy(FP_Gill[i,:], FP_Add[i,:])   
    KDE_KLD_Add[i] = stats.entropy(KDE_Gill[i,:], KDE_Add[i,:])   
  #  KLD_error_Add[i] = get_error(KDE_Gill[i,:], KDE_Add[i,:], error_Gill, error_Add[i,:], num_KDE_pts)
    
for i in range(0, num_orders):
    FP_KLD_Mult[i] = stats.entropy(FP_Gill[i,:], FP_Mult[i,:])
    KDE_KLD_Mult[i] = stats.entropy(KDE_Gill[i,:], KDE_Mult[i,:])
  #  KLD_error_Mult[i] = get_error(KDE_Gill[i,:], KDE_Mult[i,:], error_Gill, error_Mult[i,:], num_KDE_pts)
"""

"""
plt.plot(N, FP_KLD_Add_op)

plt.title("Add. noise KLD, FP")
plt.xlabel("$N$",fontsize=14)
plt.ylabel("KLD", fontsize=14)
plt.savefig("N_FP_KLD_Add.png")
plt.show()

plt.clt()

plt.plot(N, FP_KLD_Mult_op)

plt.title("Mult. noise KLD, FP")
plt.xlabel("$N$",fontsize=14)
plt.ylabel("KLD", fontsize=14)
plt.savefig("N_FP_KLD_Mult.png")
plt.show()



4


plt.plot(N, KDE_KLD_Add_op)
#plt.scatter(add_noise, KDE_KLD_Add,color='black',zorder=3)
#plt.errorbar(add_noise, KDE_KLD_Add, yerr=KLD_error_Add,fmt='none',linewidth=2,zorder=2,color='black')

plt.title("Add. noise KLD, KDE")
plt.xlabel("$N$",fontsize=14)
plt.ylabel("KLD", fontsize=14)
plt.savefig("N_KDE_KLD_Add.png")
plt.show()

plt.clt()

plt.plot(N, KDE_KLD_Mult_op)
#plt.scatter(mult_noise, KDE_KLD_Mult,color='black',zorder=3)
#plt.errorbar(mult_noise, KDE_KLD_Mult, yerr=KLD_error_Mult,fmt='none',linewidth=2,zorder=2,color='black')

plt.title("Mult. noise KLD, KDE")
plt.xlabel("$N$",fontsize=14)
plt.ylabel("KLD", fontsize=14)
plt.savefig("N_KDE_KLD_Mult.png")
plt.show()

plt.clt()
"""

P = np.zeros(tsize)
sig_a = 0
sig_m = 0
P[0] = protein_IC2

for i in range(1,tsize):
    P[i] = RK4(P[i-1], f, step_size, params)


sig_a = math.sqrt(2*dp*protein_IC)
sig_m = math.sqrt(2*dp/protein_IC)*P
KDE_Gill = bruteforce_solve(f, noise, protein_IC2, num_KDE_pts, -1, -1, params)
G1 = KDE_Gill[0]
KDE_Add = bruteforce_solve(f, g_add, protein_IC2, num_KDE_pts, -1, -1, params)
A = KDE_Add[0]
KDE_Mult = bruteforce_solve(f, g_mult, protein_IC2, num_KDE_pts, -1, -1, params)
M = KDE_Mult[0]

t = np.linspace(int_min,total_time,tsize)



#####plot noise value
plt.plot(t,A[1:])
plt.ylim(bottom=100,top=700)
plt.title("Additive noise")
plt.xlabel("time",fontsize=14)
plt.ylabel("protein", fontsize=14)
plt.savefig("add_N=1.png")
plt.show()



plt.plot(t,M[1:])
plt.ylim(bottom=100,top=700)
plt.title("Multiplicative noise")
plt.xlabel("time",fontsize=14)
plt.ylabel("protein", fontsize=14)
plt.savefig("multi_N=1.png")
plt.show()



plt.plot(t,G1[1:])
plt.ylim(bottom=100,top=700)
plt.title("Gillespie noise")
plt.xlabel("time",fontsize=14)
plt.ylabel("protein", fontsize=14)
plt.savefig("Gillespie_N=1.png")
plt.show()


plt.plot(t,P)
plt.ylim(bottom=100,top=700)
plt.title("ODE")
plt.xlabel("time",fontsize=14)
plt.ylabel("protein", fontsize=14)
plt.savefig("ODE_N=1.png")
plt.show()

    
    
#plot noise graph
p=np.linspace(0,1000,1000,True)
noiseA=np.zeros(len(p))
noiseM=np.zeros(len(p))
noiseG=np.zeros(len(p))
noise=np.zeros(len(p))

for i in range(len(p)):
    noiseA[i] = math.sqrt(2*dp*protein_IC)
    noiseM[i] = math.sqrt(2*dp/protein_IC)*p[i]
    noiseG[i] = math.sqrt((kp*kmin*G/dm)+dp*p[i])
    noise[i] = math.sqrt(2*dp*protein_IC)*((3/4)+(1/(4*protein_IC))*p[i])
    

plt.plot(p,noiseA,color='black')
plt.plot(p,noiseM,color='blue')
plt.plot(p,noiseG,color='red')
plt.plot(p,noise,color='orange')
plt.title("Noise v.s. Protein")
plt.xlabel("protein",fontsize=14)
plt.ylabel("noise", fontsize=14)
plt.savefig("noise_N=0.png")
plt.show()



# ====================================================================

# Timer
print("--- %s seconds ---" % (ti.time() - start_time))