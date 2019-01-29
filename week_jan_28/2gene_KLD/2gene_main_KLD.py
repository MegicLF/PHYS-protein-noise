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

#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt

from solution_options import num_KDE_pts
from eqn_functions import f1, f2, g1_gill, g2_gill, g1_add, g2_add, g1_mult, g2_mult
from IC_generator import generate_IC
from eqn_parameters import c


from bruteforce_solver_2D import bruteforce_solve_2D

#import sys

# =====================================================================

# Calculate +/x noise distributions and KL divergences or not. Turn these off to speed up code if just want P_gill.
do_add_comp = 1 
do_mult_comp = 1

# ============================================================================================

# siga_x, siga_y, sigm_x, sigm_y
params = [5, 5, 0.01, 0.01, c]


f = [f1, f2]
g_gill = [g1_gill, g2_gill]
g_add = [g1_add, g2_add]
g_mult = [g1_mult, g2_mult]

find_IC_step_size = 0.01
tol = 0.01
stop_tol = 0.00001

num_attractors, protein_IC = generate_IC(tol, f, find_IC_step_size, stop_tol, params)
KDE_Gill, KDE_Gill_flat, p_min, p_max = bruteforce_solve_2D(f, g_mult, protein_IC, num_attractors, num_KDE_pts, [-1,-1], [-1,-1], params)

xmin = p_min[0]
ymin = p_min[1]
xmax = p_max[0]
ymax = p_max[1]



# =======================================================================

xmax = 110
ymax = 110

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(np.rot90(KDE_Gill), cmap='viridis',  extent=[0, xmax, 0, ymax])
       
ax.set_title("$P_{ss}(x,y)$" ,fontsize=13,fontweight='normal')
ax.set_xlabel('x expression',fontsize=13)
ax.set_ylabel('y expression',fontsize=13)
ax.set_xlim([0, xmax])
ax.set_ylim([0, ymax])

basename = "mult.png"
name = "results/pss_" + basename
plt.savefig(name)
plt.show()


# ====================================================================



##############################

# ADDITIVE NOISE STUFF

if do_add_comp == 1:
    num_addnoise = 10
    add_min = 20
    add_max = 100
    add_noise = np.linspace(add_min, add_max, num_addnoise)
    
    
    KDE_Add = np.zeros((num_addnoise, num_KDE_pts, num_KDE_pts))
    KDE_Add_flat = np.zeros((num_addnoise, num_KDE_pts**2))
    
    for i in range(0, num_addnoise):
        params[0] = add_noise[i]
        KDE_Add[i,:,:], KDE_Add_flat[i,:], trash1, trash2 = bruteforce_solve_2D(f, g_add, protein_IC, num_attractors, num_KDE_pts, p_min, p_max, params)
        
    KDE_KLD_Add = np.zeros(num_addnoise)
    
    for i in range(0, num_addnoise):
        KDE_KLD_Add[i] = stats.entropy(KDE_Gill_flat, KDE_Add_flat[i,:])   
    
    plt.plot(add_noise, KDE_KLD_Add)
    plt.title("Add. noise KLD, KDE")
    plt.xlabel("$\sigma_a$",fontsize=14)
    plt.ylabel("KLD", fontsize=14)
    plt.savefig("KDE_KLD_Add.png")
    plt.show()\
    
##########################################################################
    
# MULTIPLICATIVE NOISE STUFF

if do_mult_comp == 1:
    num_multnoise = 10
    mult_min = 0.01
    mult_max = 0.1
    mult_noise = np.linspace(mult_min, mult_max, num_multnoise)
    
    KDE_Mult = np.zeros((num_multnoise, num_KDE_pts, num_KDE_pts))
    KDE_Mult_flat = np.zeros((num_multnoise, num_KDE_pts**2))
    
    for i in range(0, num_multnoise):
        params[2] = mult_noise[i]
        KDE_Mult[i,:,:], KDE_Mult_flat[i,:], trash1, trash2 = bruteforce_solve_2D(f, g_mult, protein_IC, num_attractors, num_KDE_pts, p_min, p_max, params)

    KDE_KLD_Mult = np.zeros(num_multnoise)
    
    for i in range(0, num_multnoise):
        KDE_KLD_Mult[i] = stats.entropy(KDE_Gill_flat, KDE_Mult_flat[i,:])
    
    plt.plot(mult_noise, KDE_KLD_Mult)
    plt.title("Mult. noise KLD, KDE")
    plt.xlabel("$\sigma_m$", fontsize=14)
    plt.ylabel("KLD", fontsize=14)
    plt.savefig("KDE_KLD_Mult.png")
    plt.show()




# ====================================================================

# Timer
print("--- %s seconds ---" % (ti.time() - start_time))



