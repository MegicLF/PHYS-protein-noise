# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:36:09 2018

@author: John
"""

#import matplotlib as mpl    # UNCOMMENT FOR ACCRE
#mpl.use('Agg')              # UNCOMMENT FOR ACCRE

# Timer
import time as ti
start_time = ti.time()

import numpy as np
import matplotlib.pyplot as plt

from eqn_functions import f, g_gill, hill_deriv_1
from FP_solver import FP_solve, FP_Add_solve, FP_Mult_solve
from comp_calc import KLD_calc, L2_calc, KS_calc
from rec import savegraph
import math
from eqn_parameters_og import dp

# =====================================================================

# TURN ON/OFF VARIOUS PARTS OF CODE TO SAVE TIME


do_gill_FP = 1      # Calculate Gillespie noise P_ss using FP method? 
                    # Can set to 0 if just want to recalc KLD for different range.
                 
do_add_comp = 1     # Calculuate + noise P_ss?
do_mult_comp = 0    # Calculuate x noise P_ss?

comp_type = "KLD"   # options: "KLD", "L2", "KS"       


# ============================================================================================

  
cut_off = 0.0000000001         # probability values smaller than this set to 0 in KLD calculation
naive_step = 0.1               # step size to use for Fokker-Planck integration
num_FP_pts = 200               # approximate continuous prob distribution as this many points

p_min = 1                      # minimum value to calc prob dist for; needs to be greater than 0
p_max = 1200                     # maximum value to calc prob dist for; needs to be greater than p_min
int_min = p_min/2              # lower bound of integral in Fokker-Planck calculuation;
                               # this needs to be at least several step sizes smaller than p_min

vals_FP = np.linspace(p_min, p_max, num_FP_pts)
dx = vals_FP[1] - vals_FP[0]

if do_gill_FP == 1:
    print("Calculating Gillespie noise P_ss using FP...")
    FP_Gill = FP_solve(f, g_gill, num_FP_pts, p_min, p_max, naive_step, int_min)
    #FP_Gill2 = FP_solve_OLD(f, g_gill, num_FP_pts, p_min, p_max, naive_step, int_min)
    print("Done.\n")
    
    
    plt.plot(vals_FP, FP_Gill, color='black', linewidth=3, linestyle='dashed')
    #plt.plot(vals_FP, FP_Gill2)
    plt.title("FP Gillespie $P_{ss}$")
    plt.xlabel("x",fontsize=14)
    plt.ylabel("$P_{ss}$",fontsize=14)
    savegraph("results/","FP_Gill")
    plt.show()

x0 = max(FP_Gill)



##############################

# ADDITIVE NOISE STUFF

if do_add_comp == 1:
    print("Beginning additive noise calculations.")
    
    num_addnoise = 300
    add_min = 1
    add_max = 150
    add_noise = np.linspace(add_min, add_max, num_addnoise)
    sigma_a = math.sqrt(2*dp*vals_FP[list(FP_Gill).index(x0)])
    
    grid = (1/100)*math.sqrt((sigma_a**2)/(2*abs(hill_deriv_1(x0))))
    num_pts = round((p_max-p_min)/grid)
    
    FP_Add = FP_Add_solve(num_addnoise, add_min, add_max, f, num_pts, p_min, p_max, naive_step, int_min)
    FP_Gill = FP_solve(f, g_gill, num_pts, p_min, p_max, naive_step, int_min)
    vals_FP = np.linspace(p_min, p_max, num_pts)

    plt.plot(vals_FP, FP_Gill)
    for i in range(0, num_addnoise):
        plt.plot(vals_FP, FP_Add[i,:])
    plt.show()
    
    if comp_type=="KLD":
        FP_KLD_Add = np.zeros(num_addnoise)
        for i in range(0, num_addnoise):
            FP_KLD_Add[i] = KLD_calc(FP_Gill, FP_Add[i,:], num_pts, cut_off)  
        
        x = list(FP_KLD_Add).index(min(j for j in FP_KLD_Add if j>0))
        y = add_noise[x]
        name = "("+str(round(y,5))+", "+str(round(FP_KLD_Add[x],5))+")"
        plt.plot(y,FP_KLD_Add[x],'ro')
        plt.text(y,FP_KLD_Add[x]+0.001,name)        
        plt.plot(add_noise, FP_KLD_Add)
        plt.title("Add. noise KLD, FP")
        plt.xlabel("$\sigma_a$",fontsize=14)
        plt.ylabel("KLD", fontsize=14)
        savegraph("results/","FP_KLD_Add")
        plt.show()       
    elif comp_type=="L2":
        FP_L2_Add = np.zeros(num_addnoise)
        for i in range(0, num_addnoise):
            FP_L2_Add[i] = L2_calc(FP_Gill, FP_Add[i,:], num_pts, dx)  
        
        x = list(FP_L2_Add).index(min(j for j in FP_L2_Add if j>0))
        y = add_noise[x]
        name = "("+str(round(y,5))+", "+str(round(FP_L2_Add[x],5))+")"
        plt.plot(y,FP_L2_Add[x],'ro')
        plt.text(y,FP_L2_Add[x]+0.1,name)
        plt.plot(add_noise, FP_L2_Add)
        plt.title("Add. noise L2, FP")
        plt.xlabel("$\sigma_a$",fontsize=14)
        plt.ylabel("L2", fontsize=14)
        savegraph("results/","FP_L2_Add")
        plt.show()  
    elif comp_type=="KS":
        FP_KS_Add = np.zeros(num_addnoise)
        for i in range(0, num_addnoise):
            FP_KS_Add[i] = KS_calc(FP_Gill, FP_Add[i,:], num_FP_pts)  
        
        x = list(FP_KS_Add).index(min(j for j in FP_KS_Add if j>0))
        y = add_noise[x]
        name = "("+str(round(y,5))+", "+str(round(FP_KS_Add[x],5))+")"
        plt.plot(y,FP_KS_Add[x],'ro')
        plt.text(y,FP_KS_Add[x]+0.1,name)
        plt.plot(add_noise, FP_KS_Add)
        plt.title("Add. noise KS, FP")
        plt.xlabel("$\sigma_a$",fontsize=14)
        plt.ylabel("KS", fontsize=14)
        savegraph("results/","FP_KS_Add")
        plt.show()      
    
    print("Finished additive noise calculations.\n")


    
##########################################################################
    
# MULTIPLICATIVE NOISE STUFF

if do_mult_comp == 1:
    print("Beginning multiplicative noise calculations.")
    
    num_multnoise = 100
    mult_min = 0.01
    mult_max = 0.8
    mult_noise = np.linspace(mult_min, mult_max, num_multnoise)
    
    FP_Mult = FP_Mult_solve(num_multnoise, mult_min, mult_max, f, num_FP_pts, p_min, p_max, naive_step, int_min)     

    plt.plot(vals_FP, FP_Gill)
    for i in range(0, num_multnoise):
        plt.plot(vals_FP, FP_Mult[i,:])
    plt.show()

    if comp_type=="KLD":
        FP_KLD_Mult = np.zeros(num_multnoise)
        for i in range(0, num_multnoise):
            FP_KLD_Mult[i] = KLD_calc(FP_Gill, FP_Mult[i,:], num_FP_pts, cut_off)  
        
        x = list(FP_KLD_Mult).index(min(j for j in FP_KLD_Mult if j>0))
        y = mult_noise[x]
        name = "("+str(round(y,5))+", "+str(round(FP_KLD_Mult[x],5))+")"
        plt.plot(y,FP_KLD_Mult[x],'ro')
        plt.text(y,FP_KLD_Mult[x]+0.01,name)
        
        plt.plot(mult_noise, FP_KLD_Mult)
        plt.title("Mult. noise KLD, FP")
        plt.xlabel("$\sigma_m$",fontsize=14)
        plt.ylabel("KLD", fontsize=14)
        savegraph("results/","FP_KLD_Mult")
        plt.show()       
    elif comp_type=="L2":
        FP_L2_Mult = np.zeros(num_multnoise)
        for i in range(0, num_multnoise):
            FP_L2_Mult[i] = L2_calc(FP_Gill, FP_Mult[i,:], num_FP_pts, dx)  
        
        x = list(FP_L2_Mult).index(min(j for j in FP_L2_Mult if j>0))
        y = mult_noise[x]
        name = "("+str(round(y,5))+", "+str(round(FP_L2_Mult[x],5))+")"
        plt.plot(y,FP_L2_Mult[x],'ro')
        plt.text(y,FP_L2_Mult[x]+0.01,name)
        
        plt.plot(mult_noise, FP_L2_Mult)
        plt.title("Mult. noise L2, FP")
        plt.xlabel("$\sigma_m$",fontsize=14)
        plt.ylabel("L2", fontsize=14)
        savegraph("results/","FP_L2_Mult")
        plt.show()  
    elif comp_type=="KS":
        FP_KS_Mult = np.zeros(num_multnoise)
        for i in range(0, num_multnoise):
            FP_KS_Mult[i] = KS_calc(FP_Gill, FP_Mult[i,:], num_FP_pts)  
         
        x = list(FP_KS_Mult).index(min(j for j in FP_KS_Mult if j>0))
        y = mult_noise[x]
        name = "("+str(round(y,5))+", "+str(round(FP_KS_Mult[x],5))+")"
        plt.plot(y,FP_KS_Mult[x],'ro')
        plt.text(y,FP_KS_Mult[x]+0.01,name)
        
        plt.plot(mult_noise, FP_KS_Mult)
        plt.title("Mult. noise KS, FP")
        plt.xlabel("$\sigma_m$",fontsize=14)
        plt.ylabel("KS", fontsize=14)
        savegraph("results/","FP_KS_Mult")
        plt.show()  
    
    print("Finished multiplicative noise calculations.\n")


# ====================================================================

# Timer
print("--- %s seconds ---" % (ti.time() - start_time))



