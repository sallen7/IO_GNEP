######### 1_generating_more_graphs_different_costs_sioux_falls.py #################

############ Script to Create the Box Plots #############

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
import csv

############ Creating the Box Plots ##################
num_trials = 10
type_of_alpha = [1,2]
number_of_players = [2,5,10]

#############################################
flow_error = []
min_flow_error = []
median_flow_error = []
max_flow_error = []

flow_error_normalized = []
min_flow_error_normalized = []
median_flow_error_normalized = []
max_flow_error_normalized = []

timing = []
min_timing = []
median_timing = []
max_timing = []

objective_error = []
min_objective_error = []
median_objective_error = []
max_objective_error = []

names = []
######## Aggregating the Data ##############
for player in number_of_players:
    for alpha in type_of_alpha:
        if alpha == 1:
            alpha_var = int(round(player/2))
            alpha_graph = player/2
        elif alpha == 2:
            alpha_var = int(player)
            alpha_graph = player
        
        if player == 10 and alpha == 1:
            pass
        else:        
            ############ Transitioning to Sioux Falls ##############
            name_of_data = "data_saving_experiment_different_costs_Sioux_Falls"+\
                "_num_players_"+str(player)+"_alpha_"+str(alpha_var)+".csv"
            names.append(str(player)+"/"+str(alpha_graph))
            data = pandas.read_csv(name_of_data)
                        
            #pdb.set_trace()
            
            ph_flow_error = np.zeros((num_trials,))
            ph_flow_error_normalized = np.zeros((num_trials,))
            ph_timing = np.zeros((num_trials,))
            ph_obj_error = np.zeros((num_trials,))
            normalization_factor = (player)*(76)*(552)
            
            #pdb.set_trace()
            for j in range(0,num_trials):
                ph_flow_error[j] = float(data.loc[2,str(j+1)])
                ph_flow_error_normalized[j] = float(data.loc[2,str(j+1)])/(normalization_factor)
                ph_timing[j] = (float(data.loc[3,str(j+1)])/60)
                ph_obj_error[j] = float(data.loc[4,str(j+1)])
                if data.loc[5,str(j+1)]=="optimal":
                    pass
                else:
                    print("PROBLEM: Non-optimal flag")
                    pdb.set_trace()
            
            ########## Calculating the Min, Median, and Max ############
            min_flow_error.append(np.min(ph_flow_error))
            median_flow_error.append(np.median(ph_flow_error))
            max_flow_error.append(np.max(ph_flow_error))
            
            min_flow_error_normalized.append(np.min(ph_flow_error_normalized))
            median_flow_error_normalized.append(np.median(ph_flow_error_normalized))
            max_flow_error_normalized.append(np.max(ph_flow_error_normalized))
            
            min_timing.append(np.min(ph_timing))
            median_timing.append(np.median(ph_timing))
            max_timing.append(np.max(ph_timing))
        
            min_objective_error.append(np.min(ph_obj_error))
            median_objective_error.append(np.median(ph_obj_error))
            max_objective_error.append(np.max(ph_obj_error))    
            
            ######## Appending to the Lists ##############
            flow_error.append(ph_flow_error)
            flow_error_normalized.append(ph_flow_error_normalized)
            timing.append(ph_timing)
            objective_error.append(ph_obj_error)

########### Creating the Box Plots #############
##### Print out Min/Median/Max ######
print("min flow error", min_flow_error)
print("median_flow_error",median_flow_error)
print("max_flow_error",max_flow_error)
print("*********************************")
print("min flow error normalized", min_flow_error_normalized)
print("median_flow_error normalized",median_flow_error_normalized)
print("max_flow_error normalized",max_flow_error_normalized)
print("******************************")
print("min_timing",min_timing)
print("median_timing",median_timing)
print("max_timing",max_timing)
print("*****************************")
print("min_objective_error",min_objective_error)
print("median_objective_error",median_objective_error)
print("max_objective_error",max_objective_error)

##### Flow Error #########
visual,test = plt.subplots()
test.boxplot(flow_error)
test.set_title("Flow Error for SF: Different Parameters and 10 Trials for Each Experiment")
#pdb.set_trace()
plt.xticks(rotation=45)
test.set_xticklabels(names)
test.set_ylabel("Frobenius-norm")
    
plt.show()


##### Flow Error Normalized #########
test1 = plt.subplot(111)
test1.boxplot(flow_error_normalized)
plt.title("Normalized Flow Error for SF: Different Parameters and 10 Trials for Each Experiment")
plt.xticks(rotation=45)
test1.set_xticklabels(names)
plt.ylabel("Frobenius-norm/((OD Pair Number)*(Num Arcs)*(Num Players))")
plt.show()


##### Timing ##########
test2 = plt.subplot(111)
test2.boxplot(timing)
plt.title("Timing for SF: Different Parameters and 10 Trials for Each Experiment")
plt.xticks(rotation=45)
test2.set_xticklabels(names)
plt.ylabel("Minutes")
plt.show()

#### Objective Error ######
test3 = plt.subplot(111)
test3.boxplot(objective_error)
plt.title("Objective Function Values for SF: Different Parameters and 10 Trials for Each Experiment")
plt.xticks(rotation=45)
test3.set_xticklabels(names)
plt.ylabel("Objective Function Value")
plt.show()

