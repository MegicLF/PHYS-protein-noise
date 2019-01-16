# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:26:15 2018

@author: John
"""



# Set up simulation.


# Fokker-Planck solution options.
#-----------------------------------
naive_step = 0.1
int_min = 10
num_FP_pts = 100






# SDE simulation / KDE options.
#-----------------------------------
tsize = 300                                     # number of time steps
total_time = 25                                 # total simulation time   
step_size = total_time/tsize                    # euler-maruyama step size
 
perc_use = 0.9                                # Percent of time steps to use from latter part of each simulation
start_index = int(tsize - perc_use*tsize)       # Index corresponding to the specified percent of time steps to use

num_sims = 20                                   # Number of simulations  

num_KDE_pts = 100






#tsize = 200                                     # number of time steps
#total_time = 10                                 # total simulation time   
#step_size = total_time/tsize                    # euler-maruyama step size
# 
#perc_use = 0.6                                  # Percent of time steps to use from latter part of each simulation
#start_index = int(tsize - perc_use*tsize)       # Index corresponding to the specified percent of time steps to use
#
#num_sims = 40                                   # Number of simulations  
#
#num_KDE_pts = 100