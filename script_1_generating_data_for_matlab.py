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

################### Step 1: Generating Data ###############3#######

######################### nxn Grid ####################################

nxn_grid = nx.generators.lattice.grid_2d_graph(int(sys.argv[1]),int(sys.argv[1]))
incidence_matrix = nx.linalg.graphmatrix.incidence_matrix(nxn_grid)
incid_mat = incidence_matrix.todense()
(num_nodes,num_arcs) = np.shape(incid_mat)
for i in range(0,num_arcs):
    ph = incid_mat[:,i]
    #for j in range(0,num_nodes):
    j = 0
    while ph[j] != 1:
        j = j + 1
    incid_mat[j,i] = -1

full_incidence_matrix = np.concatenate((incid_mat,-1*incid_mat),1)

###################################################################################
################### Step 2: Setting up Object and Saving Matlab #############################

name_of_grid = str(sys.argv[1])+"x"+str(sys.argv[1])+"_Grid"

GRKF_Object = generalized_RK_framework(num_nodes=num_nodes,num_arcs=num_arcs*2,num_players=int(sys.argv[3]),num_trials=10,\
                                       node_arc_incidence_matrix=full_incidence_matrix,\
                                       name_of_graph=name_of_grid)

alpha_flag = int(sys.argv[2])
if alpha_flag == 1:
    alpha = float(sys.argv[3])*0.5
elif alpha_flag == 2:
    alpha = float(sys.argv[3])
    
    
GRKF_Object.saving_for_matlab_files_randomized_costs(lowerbound_c=1,upperbound_c=5,\
                                                 lowerbound_chat=5,upperbound_chat=20,\
                                                 alpha=alpha,if_different_costs=1)

################### Step 3: Saving the Object #################################
#https://www.datacamp.com/community/tutorials/pickle-python-tutorial

name_of_file = "class_object_1"
test = open(name_of_file,'wb')
pickle.dump(GRKF_Object,test)
test.close()

#https://www.mathworks.com/matlabcentral/answers/327116-run-function-from-command-line
#https://www.mathworks.com/matlabcentral/answers/410079-how-can-i-execute-an-m-file-from-windows-command-line-without-opening-matlab-s-command-window
#https://www.mathworks.com/matlabcentral/answers/479672-how-to-run-batch-file-in-matlab
#^The site that helped with the MATLAB command line code 
