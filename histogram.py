# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 09:35:59 2018

@author: fengl
"""

import parameter
import numpy as np
import matplotlib.pyplot as plt
import QSS
import steadystate
import pandas as pd

steadystate

initial_p1 = steadystate.my_ss2
initial_p2 = steadystate.my_ss2
    

N = parameter.Num_sim  #number of simulations

p_A1 = [[0] for a in range(N)]
p_M1 = [[0] for a in range(N)]
p_G1 = [[0] for a in range(N)]
p_A2 = [[0] for a in range(N)]
p_M2 = [[0] for a in range(N)]
p_G2 = [[0] for a in range(N)]
p_A = [[0] for a in range(N)]
p_M = [[0] for a in range(N)]
p_G = [[0] for a in range(N)]


def write(initial_p,i,j):
    name = 'A='+str(round(parameter.p_AddNoise[j],2))+'M='+str(round(parameter.p_MultiNoise[j],2))+','+str(round(initial_p2,3))+" ("+str(i)+").csv"
    dataframe = pd.DataFrame({'protein(Gillespie noise)':QSS.p_Gillespie,
                              'protein(multiplicative noise)':QSS.p_multi,
                              'protein(additive noise)':QSS.p_add,  
                              'time':QSS.t})
    dataframe.to_csv(name,index=False,sep=',')
   
for j in range (parameter.noise_size):    
    if abs(initial_p1-initial_p2) > 0.1:
        for i in range(N):
            QSS.QSS(initial_p1,parameter.p_AddNoise[j],parameter.p_MultiNoise[j])
            write(initial_p1,i,j)
            p_A1[i] = QSS.p_add[parameter.data_size:]
            p_M1[i] = QSS.p_multi[parameter.data_size:]
            p_G1[i] = QSS.p_Gillespie[parameter.data_size:]
    
        for i in range(N):
            QSS.QSS(initial_p2,parameter.p_AddNoise[j],parameter.p_MultiNoise[j])
            write(initial_p2,i,j)
            p_A2[i] = QSS.p_add[parameter.data_size:]
            p_M2[i] = QSS.p_multi[parameter.data_size:]
            p_G2[i] = QSS.p_Gillespie[parameter.data_size:]
        
        Add1 = np.array(p_A1[0])
        Multi1 = np.array(p_M1[0])
        Gille1 = np.array(p_G1[0])
        Add2 = np.array(p_A2[0])
        Multi2 = np.array(p_M2[0])
        Gille2 = np.array(p_G2[0])
        
        for k in range(N):
            Add1 = np.append(Add1,p_A1[k])
            Multi1 = np.append(Multi1,p_M1[k])
            Gille1 = np.append(Gille1,p_G1[k])
            
        for k in range(N):
            Add2 = np.append(Add2,p_A2[k])
            Multi2 = np.append(Multi2,p_M2[k])
            Gille2 = np.append(Gille2,p_G2[k])
            

        t1 = 'A'+str(round(parameter.p_AddNoise[j],2))+','+str(round(initial_p1,2))+'.png'
        plt.hist(Add1)
        plt.title("additive noise: "+str(round(initial_p1,2)))
        plt.savefig(t1)
        plt.show()
        plt.hist(Multi1)
        t2 = 'M'+str(round(parameter.p_MultiNoise[j],2))+','+str(round(initial_p1,2))+'.png'
        plt.title("multiple noise: "+str(round(initial_p1,2)))
        plt.savefig(t2)
        plt.show()
        plt.hist(Gille1)
        t3 = 'G'+str(round(initial_p1,2))+'.png'
        plt.title("Gillespie noise: "+str(round(initial_p1,2)))
        plt.savefig(t3)
        plt.show()
        
        t4 = 'A'+str(round(parameter.p_AddNoise[j],2))+','+str(round(initial_p2,2))+'.png'
        plt.hist(Add2)
        plt.title("additive noise: "+str(round(initial_p2,2)))
        plt.savefig(t4)
        plt.show()
        plt.hist(Multi2)
        t5 = 'M'+str(round(parameter.p_MultiNoise[j],2))+','+str(round(initial_p2,2))+'.png'
        plt.title("multiple noise: "+str(round(initial_p2,2)))
        plt.savefig(t5)
        plt.show()
        plt.hist(Gille2)
        t6 = 'G'+str(round(initial_p2,2))+'.png'
        plt.title("Gillespie noise: "+str(round(initial_p2,2)))
        plt.savefig(t6)
        plt.show()

        
    else:
        initial_p = (initial_p1+initial_p2)/2
        for i in range(N):
            QSS.QSS(initial_p,parameter.p_AddNoise[j],parameter.p_MultiNoise[j])
            write(initial_p,i,j)
            p_A[i] = QSS.p_add[parameter.data_size:]
            p_M[i] = QSS.p_multi[parameter.data_size:]
            p_G[i] = QSS.p_Gillespie[parameter.data_size:]
        
        Add = np.array(p_A[0])
        Multi = np.array(p_M[0])
        Gille = np.array(p_G[0])
    
        
        for k in range(N):
            Add = np.append(Add,p_A[k])
            Multi = np.append(Multi,p_M[k])
            Gille = np.append(Gille,p_G[k])

        t1 = 'A'+str(round(parameter.p_AddNoise[j],2))+','+str(round(initial_p,2))+'.png'
        plt.hist(Add)
        plt.title("additive noise: "+str(round(initial_p,2)))
        plt.savefig(t1)
        plt.show()
        plt.hist(Multi)
        t2 = 'M'+str(round(parameter.p_MultiNoise[j],2))+','+str(round(initial_p,2))+'.png'
        plt.title("multiple noise: "+str(round(initial_p,2)))
        plt.savefig(t2)
        plt.show()
        plt.hist(Gille)
        t3 = 'G'+str(round(initial_p,2))+'.png'
        plt.title("Gillespie noise: "+str(round(initial_p,2)))
        plt.savefig(t3)
        plt.show()

