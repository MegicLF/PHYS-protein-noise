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

import matplotlib as mpl           # uncomment these to run on ACCRE
mpl.use('Agg')
import matplotlib.pyplot as plt

from solution_options import naive_step, int_min, num_KDE_pts, num_FP_pts
from eqn_functions import f, noise, noise_deriv, g_add, g_prime_add, g_mult, g_prime_mult
from get_IC import get_IC 
from eqn_parameters import kp, G, kmin, kmax, dp, dm, Nmin, Nmax


from FP_solver import FP_solve
from bruteforce_solver import bruteforce_solve

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

num_sigmas = 100
sig_a_min = 1
sig_a_max = 100
sig_a = np.linspace(sig_a_min, sig_a_max, num_sigmas)

sig_m_min = 0.01
sig_m_max = 2
sig_m = np.linspace(sig_m_min, sig_m_max, num_sigmas)
epsilon = 0



# IC lower/upper bounds and RK4 step size. 
pss_lowerbound = (kp*kmin*G)/(dp*dm)            # lower bound for p_ss
pss_upperbound = (kp*kmax*G)/(dp*dm)            # upper bound for p_ss
get_IC_stepsize = 0.01                          # step size used when finding IC  


KDE_KLDmin_a = 1000*np.ones(num_orders)
KDE_KLDmin_m = 1000*np.ones(num_orders)
FP_KLDmin_a = 1000*np.ones(num_orders)
FP_KLDmin_m = 1000*np.ones(num_orders)

KDE_best_sig_a = np.zeros(num_orders)
KDE_best_sig_m = np.zeros(num_orders)
FP_best_sig_a = np.zeros(num_orders)
FP_best_sig_m = np.zeros(num_orders)

for j in range(0, num_orders):
    # 1: sig_a, 2: sig_m, 3: epsilon, 4: order of feedback N
    params = [sig_a[0], sig_m[0], epsilon, N[j]]


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
    else:
        # Probably 2 steady states.
        num_steadystates = 2    
        protein_IC = [pss_low, pss_high]
        
    # -------------------------------------------------------------------------------------------------    
    
    
    
    KDE_Gill[j,:], p_min, p_max = bruteforce_solve(f, noise, protein_IC, num_KDE_pts, -1, -1, params)
    #FP_Gill[j,:] = FP_solve(f, noise, noise_deriv, num_FP_pts, p_min, p_max, naive_step, int_min, params)
    
    
    for i in range(0, num_sigmas):
        params = [sig_a[i], sig_m[i], epsilon, N[j]]
        
        
        KDE_Add[j,:], trash1, trash2 = bruteforce_solve(f, g_add, protein_IC, num_KDE_pts, p_min, p_max, params)
        #FP_Add[j,:] = FP_solve(f, g_add, g_prime_add, num_FP_pts, p_min, p_max, naive_step, int_min, params)
    
        KDE_Mult[j,:], trash1, trash2 = bruteforce_solve(f, g_mult, protein_IC, num_KDE_pts, p_min, p_max, params)
        #FP_Mult[j,:] = FP_solve(f, g_mult, g_prime_mult, num_FP_pts, p_min, p_max, naive_step, int_min, params)
        
        
        
        KDE_KLDtry_a = stats.entropy(KDE_Gill[j,:], KDE_Add[j,:])   
        print(KDE_KLDtry_a)
        #FP_KLDtry_a = stats.entropy(FP_Gill[j,:], FP_Add[j,:])  
        
        KDE_KLDtry_m = stats.entropy(KDE_Gill[j,:], KDE_Mult[j,:])   
        #FP_KLDtry_m = stats.entropy(FP_Gill[j,:], FP_Mult[j,:])   
        
        if KDE_KLDtry_a < KDE_KLDmin_a[j]:
            KDE_KLDmin_a[j] = KDE_KLDtry_a
            KDE_best_sig_a[j] = sig_a[i]
        
        #if FP_KLDtry_a > FP_KLDmin_a[j]:
         #   FP_KLDmin_a[j] = FP_KLDtry_a
            
        if KDE_KLDtry_m < KDE_KLDmin_m[j]:
            KDE_KLDmin_m[j] = KDE_KLDtry_m
            KDE_best_sig_m[j] = sig_m[i]
            
        #if FP_KLDtry_m > KDE_KLDmin_m[j]:
         #   FP_KLDmin_m[j] = FP_KLDtry_m
        
    ###
    #test = 1
# ============================================================================================





vals_KDE = np.linspace(p_min, p_max, num_KDE_pts)
vals_FP = np.linspace(p_min, p_max, num_FP_pts)


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


#FP_KLD_Add = np.zeros(num_orders)
#FP_KLD_Mult = np.zeros(num_orders)
#
#KDE_KLD_Add = np.zeros(num_orders)
#KDE_KLD_Mult = np.zeros(num_orders)
#
#KLD_error_Add = np.zeros(num_orders)
#KLD_error_Mult = np.zeros(num_orders)
#
#
#for i in range(0, num_orders):
#    FP_KLD_Add[i] = stats.entropy(FP_Gill[i,:], FP_Add[i,:])   
#    KDE_KLD_Add[i] = stats.entropy(KDE_Gill[i,:], KDE_Add[i,:])   
#  #  KLD_error_Add[i] = get_error(KDE_Gill[i,:], KDE_Add[i,:], error_Gill, error_Add[i,:], num_KDE_pts)
#    
#for i in range(0, num_orders):
#    FP_KLD_Mult[i] = stats.entropy(FP_Gill[i,:], FP_Mult[i,:])
#    KDE_KLD_Mult[i] = stats.entropy(KDE_Gill[i,:], KDE_Mult[i,:])
#  #  KLD_error_Mult[i] = get_error(KDE_Gill[i,:], KDE_Mult[i,:], error_Gill, error_Mult[i,:], num_KDE_pts)



plt.plot(N, KDE_KLDmin_a)
plt.title("Best KLD (+ noise) vs N, KDE")
plt.xlabel("$N$",fontsize=14)
plt.ylabel("Best KLD", fontsize=14)
plt.savefig("results/N_KDE_KLD_Add.png")
plt.show()
plt.gcf().clear()

plt.plot(N, KDE_best_sig_a)
plt.title("Best $\sigma$ (+ noise) vs N, KDE")
plt.xlabel("$N$",fontsize=14)
plt.ylabel("Best $\sigma$", fontsize=14)
plt.savefig("results/N_KDE_bestsig_Add.png")
plt.show()
plt.gcf().clear()


#plt.plot(N, FP_KLDmin_a)
#plt.title("Best KLD (+ noise) vs N, FP")
#plt.xlabel("$N$",fontsize=14)
#plt.ylabel("Best KLD", fontsize=14)
#plt.savefig("results/N_FP_KLD_Add.png")
#plt.show()
#plt.gcf().clear()


plt.plot(N, KDE_KLDmin_m)
plt.title("Best KLD (x noise) vs N, KDE")
plt.xlabel("$N$",fontsize=14)
plt.ylabel("Best KLD", fontsize=14)
plt.savefig("results/N_KDE_KLD_Mult.png")
plt.show()
plt.gcf().clear()

plt.plot(N, KDE_best_sig_m)
plt.title("Best $\sigma$ (x noise) vs N, KDE")
plt.xlabel("$N$",fontsize=14)
plt.ylabel("Best $\sigma$", fontsize=14)
plt.savefig("results/N_KDE_bestsig_Mult.png")
plt.show()
plt.gcf().clear()


#plt.plot(N, FP_KLDmin_m)
#plt.title("Best KLD (x noise) vs N, FP")
#plt.xlabel("$N$",fontsize=14)
#plt.ylabel("Best KLD", fontsize=14)
#plt.savefig("results/N_FP_KLD_Mult.png")
#plt.show()
#plt.gcf().clear()







# ====================================================================

# Timer
print("--- %s seconds ---" % (ti.time() - start_time))



