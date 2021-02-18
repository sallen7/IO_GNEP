########## Script 1 ###################

import sys

from RK_IO_model import RK_IO_methods
from Generalized_RK_Framework import generalized_RK_framework

import pdb #for debugging
import numpy as np                     
import pyomo.environ as pyo
from pyomo.opt import SolverFactory 
from pyomo.opt import SolverStatus, TerminationCondition
import pyomo.mpec as pyompec #for the complementarity
import math
from scipy.io import savemat, loadmat
import pandas
import time
import matplotlib.pyplot as plt
import pickle
import networkx as nx

################### Step 1: Generating Data ######################

##################### Sioux Falls ###################################
##### Thanks to the Github Jupyter code that came with this
##### and the pandas documentation (plus this site https://www.geeksforgeeks.org/indexing-and-selecting-data-with-pandas/)

sioux_falls_network = pandas.read_csv("SiouxFalls_flow.tntp",\
                                      sep="\t")

incidence_matrix = np.zeros((24,76))

for i in range(0,76):
    end = sioux_falls_network.loc[i,"To "]
    start = sioux_falls_network.loc[i,"From "]
    incidence_matrix[end-1,i] = 1
    incidence_matrix[start-1,i] = -1
    
###################################################################################
################### Step 2: Setting up Object and Saving Matlab #############################

name_of_grid = "Sioux_Falls"


GRKF_Object = generalized_RK_framework(num_nodes=24,num_arcs=76,num_players=int(sys.argv[2]),num_trials=10,\
                                       node_arc_incidence_matrix=incidence_matrix,\
                                       name_of_graph=name_of_grid)

alpha_flag = int(sys.argv[1])
if alpha_flag == 1:
    alpha = float(sys.argv[2])*0.5
elif alpha_flag == 2:
    alpha = float(sys.argv[2])
    
    
GRKF_Object.saving_for_matlab_files_randomized_costs(lowerbound_c=1,upperbound_c=5,\
                                                 lowerbound_chat=5,upperbound_chat=20,\
                                                 alpha=alpha,if_different_costs=0)

################### Step 3: Saving the Object #################################
#https://www.datacamp.com/community/tutorials/pickle-python-tutorial

name_of_file = "class_object_1"
test = open(name_of_file,'wb')
pickle.dump(GRKF_Object,test)
test.close()

