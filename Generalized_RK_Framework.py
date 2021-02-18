####### Class for Python Files in RK Framework ###########

from RK_IO_model import RK_IO_methods

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

class generalized_RK_framework():
    
    def __init__(self,num_nodes,num_arcs,num_players,num_trials,node_arc_incidence_matrix,\
                 name_of_graph):
        # The incidence matrix needs to be a numpy array
        
        self.num_nodes = num_nodes
        self.num_arcs = num_arcs
        self.num_players = num_players
        self.num_trials = num_trials
        self.incidence_matrix = node_arc_incidence_matrix
        self.name_of_graph = name_of_graph
        
    def saving_for_matlab_files_randomized_costs(self,lowerbound_c,upperbound_c,\
                                                 lowerbound_chat,upperbound_chat,\
                                                 alpha,if_different_costs):
        # This method works when you dont have any specific costs that you are passing
        
        self.lowerbound_c = lowerbound_c
        self.upperbound_c = upperbound_c
        self.lowerbound_chat = lowerbound_chat
        self.upperbound_chat = upperbound_chat
        self.alpha = alpha
        self.if_different_costs = if_different_costs
        
        savemat("data_for_matlab_original_costs.mat",{"node_arc_incidence":self.incidence_matrix,\
                                                      "num_nodes":self.num_nodes,\
                                                      "num_arcs":self.num_arcs,\
                                                      "num_players":self.num_players,\
                                                      "num_trials":self.num_trials,\
                                                      "lowerbound_c":lowerbound_c,\
                                                      "upperbound_c":upperbound_c,\
                                                      "lowerbound_chat":lowerbound_chat,\
                                                      "upperbound_chat":upperbound_chat,\
                                                      "alpha":alpha})
    
    def running_IO_code_to_obtain_costs(self):
        
        ############# Step 0: Establish the CSV Files ##############
        original_cost_flow_values = {}
        for i in range(1,self.num_trials+1):
            original_cost_flow_values[i] = "traffic_results_original_costs_"+str(i)+".csv"
        
        original_cost_mat_files = {}
        for i in range(1,self.num_trials+1):
            original_cost_mat_files[i] = "costs_iteration_"+str(i)+".mat"
        
        
        ########### Error and Timing Dictionaries ###########
        errors_c = {}
        errors_chat = {}
        timing = {}
        obj_func_values = {}
        return_values = {}
        
        
        ######### Figuring out number of OD pairs ##########
        dummy_variable = 1
        for i in range(1,self.num_nodes+1):
            dummy_variable = dummy_variable*i
        
        dummy_variable2 = 1
        for i in range(1,self.num_nodes-2+1):
            dummy_variable2 = dummy_variable2*i
            
        num_od_pairs = int((dummy_variable/(2*dummy_variable2))*2)
        
        
        ############# Step 1: Start the Trials ######################
        for iteration in range(1,self.num_trials+1):
            print("**************** Iteration ", iteration, "************************")
            
            ################ Step 2: Bring in the Flow Data ################### 
            
            counter = 1
            equilibrium_data_dict = {}
            for i in range(0,(self.num_arcs*self.num_players+1)*num_od_pairs,\
                           self.num_arcs*self.num_players+1):
                equilibrium_data = pandas.read_csv(original_cost_flow_values[iteration],\
                                                   skiprows=i,nrows=self.num_arcs*self.num_players+1)
                sets = []
                for j in range(0,self.num_arcs*self.num_players,self.num_arcs):
                    sets.append((j,j+self.num_arcs))
                
                ph = {}
                counter2 = 1
                for (j,k) in sets:
                    for g in range(j,k):
                        ph[counter2,g+1-j] = float(equilibrium_data.loc[g,'data'])
                    counter2 = counter2 + 1
                
                equilibrium_data_dict[counter] = ph
                counter = counter + 1
            
            ############ Step 3: Solving the IO/RK Model ##############
    
            begin = time.perf_counter()
            RK_instance = RK_IO_methods(data=equilibrium_data_dict,\
                            amount_of_data=num_od_pairs,number_of_nodes=self.num_nodes,\
                            number_of_arcs=self.num_arcs,incidence_matrix=self.incidence_matrix,\
                            lowerbound=self.lowerbound_c,upperbound=self.upperbound_c,\
                            lowerbound_c_hat=self.lowerbound_chat,\
                            upperbound_c_hat=self.upperbound_chat,alpha=self.alpha)
            
            RK_instance.RK_model_n_players(n=self.num_players)
            RK_instance.RK_model_n_player_solve(epsilon = 1e-12) 
            end = time.perf_counter()
        
            timing[iteration] = end - begin
            
            obj_func_values[iteration] = RK_instance.obj_func_value
            return_values[iteration] = RK_instance.return_value
            
            ############ Step 4: Calculating the Errors #############
            
            IO_c_values = RK_instance.c_values
            IO_chat_values = RK_instance.c_hat_values
            
            original_costs = loadmat(original_cost_mat_files[iteration])
            original_c = original_costs['c_vector']
            original_chat = original_costs['c_hat_vector']
            
            ph1 = 0
            ph2 = 0
            for i in range(1,self.num_arcs+1):
                ph1 = ph1 + (IO_c_values[i]-original_c[0,i-1])**2
                ph2 = ph2 + (IO_chat_values[i] - original_chat[0,i-1])**2
            
            errors_c[iteration] = math.sqrt(ph1)
            errors_chat[iteration] = math.sqrt(ph2)
            
            print("errors_c",errors_c)
            print("errors_chat",errors_chat)
            print("timing (seconds)",timing)
            print("obj_func_values",obj_func_values)
            
            ############### Step 5: Saving the Files ###########
            IO_costs_c_vec = np.zeros((self.num_arcs,))
            IO_costs_chat_vec = np.zeros((self.num_arcs,))
            
            for i in range(0,self.num_arcs):
                IO_costs_c_vec[i] = IO_c_values[i+1]
                IO_costs_chat_vec[i] = IO_chat_values[i+1]
            
            savemat("IO_costs_"+str(iteration)+".mat",\
                    {"IO_costs_c_vec":IO_costs_c_vec,\
                     "IO_costs_c_hat_vec":IO_costs_chat_vec})
            
            del RK_instance
        
        
        ############## Step 6: Saving the Errors ###########
        self.errors_c = errors_c
        self.errors_chat = errors_chat
        self.timing = timing
        self.obj_func_values = obj_func_values
        self.return_values = return_values

    def running_IO_code_to_obtain_costs_different_costs(self):
        
        ############# Step 0: Establish the CSV Files ##############
        original_cost_flow_values = {}
        for i in range(1,self.num_trials+1):
            original_cost_flow_values[i] = "traffic_results_original_costs_"+str(i)+".csv"
        
        original_cost_mat_files = {}
        for i in range(1,self.num_trials+1):
            original_cost_mat_files[i] = "costs_iteration_"+str(i)+".mat"
        
        
        ########### Error and Timing Dictionaries ###########
        errors_c = {}
        errors_chat = {}
        timing = {}
        obj_func_values = {}
        return_values = {}
        
        ######### Figuring out number of OD pairs ##########
        dummy_variable = 1
        for i in range(1,self.num_nodes+1):
            dummy_variable = dummy_variable*i
        
        dummy_variable2 = 1
        for i in range(1,self.num_nodes-2+1):
            dummy_variable2 = dummy_variable2*i
            
        num_od_pairs = int((dummy_variable/(2*dummy_variable2))*2)
        
        ############# Step 1: Start the Trials ######################
        for iteration in range(1,self.num_trials+1):
            print("**************** Iteration ", iteration, "************************")
            
            ################ Step 2: Bring in the Flow Data ################### 
            
            counter = 1
            equilibrium_data_dict = {}
            for i in range(0,(self.num_arcs*self.num_players+1)*num_od_pairs,\
                           self.num_arcs*self.num_players+1):
                equilibrium_data = pandas.read_csv(original_cost_flow_values[iteration],\
                                                   skiprows=i,nrows=self.num_arcs*self.num_players+1)
                #pdb.set_trace()
                sets = []
                for j in range(0,self.num_arcs*self.num_players,self.num_arcs):
                    sets.append((j,j+self.num_arcs))
                
                ph = {}
                counter2 = 1
                for (j,k) in sets:
                    for g in range(j,k):
                        ph[counter2,g+1-j] = float(equilibrium_data.loc[g,'data'])
                    counter2 = counter2 + 1
                
                equilibrium_data_dict[counter] = ph
                counter = counter + 1
            
            ############ Step 3: Solving the IO/RK Model ##############
    
            begin = time.perf_counter()
            RK_instance = RK_IO_methods(data=equilibrium_data_dict,\
                            amount_of_data=num_od_pairs,number_of_nodes=self.num_nodes,\
                            number_of_arcs=self.num_arcs,incidence_matrix=self.incidence_matrix,\
                            lowerbound=self.lowerbound_c,upperbound=self.upperbound_c,\
                            lowerbound_c_hat=self.lowerbound_chat,\
                            upperbound_c_hat=self.upperbound_chat,alpha=self.alpha)
            
            RK_instance.RK_model_n_players_different_costs(n=self.num_players)
            RK_instance.RK_model_n_player_different_costs_solve(epsilon = 1e-12) 
            end = time.perf_counter()
        
            timing[iteration] = end - begin
            obj_func_values[iteration] = RK_instance.obj_func_value
            return_values[iteration] = RK_instance.return_value
            
            ############ Step 4: Calculating the Errors #############
            
            IO_c_mat = RK_instance.c_matrix 
            IO_chat_mat = RK_instance.c_hat_matrix
            
            original_costs = loadmat(original_cost_mat_files[iteration])
            original_c_mat = original_costs['c_matrix']
            original_chat_mat = original_costs['c_hat_matrix']
            
            ph1 = 0
            ph2 = 0
            for p in range(1,self.num_players+1):
                for i in range(1,self.num_arcs+1):
                    ph1 = ph1 + (IO_c_mat[p,i]-original_c_mat[p-1,i-1])**2
                    ph2 = ph2 + (IO_chat_mat[p,i] - original_chat_mat[p-1,i-1])**2
            
            errors_c[iteration] = math.sqrt(ph1)
            errors_chat[iteration] = math.sqrt(ph2)
            
            print("obj_func_value",obj_func_values)
            print("errors_c",errors_c)
            print("errors_chat",errors_chat)
            print("timing (seconds)",timing)
            
            ############### Step 5: Saving the Files ###########
            IO_costs_c_mat = np.zeros((self.num_players,self.num_arcs))
            IO_costs_chat_mat = np.zeros((self.num_players,self.num_arcs))
            
            for p in range(0,self.num_players):
                for i in range(0,self.num_arcs):
                    IO_costs_c_mat[p,i] = IO_c_mat[p+1,i+1]
                    IO_costs_chat_mat[p,i] = IO_chat_mat[p+1,i+1]
            
            savemat("IO_costs_"+str(iteration)+".mat",\
                    {"IO_costs_c_mat":IO_costs_c_mat,\
                     "IO_costs_c_hat_mat":IO_costs_chat_mat})
            
            RK_instance = None
        
        
        ############## Step 6: Saving the Errors ###########
        self.errors_c = errors_c
        self.errors_chat = errors_chat
        self.timing = timing
        self.obj_func_values = obj_func_values
        self.return_values = return_values
        
    def calculating_flow_error(self):
        
        ############# Step 0: Create the File Names ###############
        original_cost_flow_values = {}
        for i in range(1,self.num_trials+1):
            original_cost_flow_values[i] = "traffic_results_original_costs_"+str(i)+".csv"
        
        IO_cost_flow_values = {}
        for i in range(1,self.num_trials+1):
            IO_cost_flow_values[i] = "IO_traffic_results_"+str(i)+".csv"
        
        differences = {}
        ######### Figuring out number of OD pairs ##########
        dummy_variable = 1
        for i in range(1,self.num_nodes+1):
            dummy_variable = dummy_variable*i
        
        dummy_variable2 = 1
        for i in range(1,self.num_nodes-2+1):
            dummy_variable2 = dummy_variable2*i
            
        num_od_pairs = int((dummy_variable/(2*dummy_variable2))*2)
        
        
        for iteration in range(1,self.num_trials+1):
            #################### Original Flow Values ########################
            counter = 1
            equilibrium_data_dict = {}
            for i in range(0,(self.num_arcs*self.num_players+1)*num_od_pairs,\
                           self.num_arcs*self.num_players+1):
                equilibrium_data = pandas.read_csv(original_cost_flow_values[iteration],\
                                                   skiprows=i,nrows=self.num_arcs*self.num_players+1)
                #pdb.set_trace()
                sets = []
                for j in range(0,self.num_arcs*self.num_players,self.num_arcs):
                    sets.append((j,j+self.num_arcs))
                
                ph = {}
                counter2 = 1
                for (j,k) in sets:
                    for g in range(j,k):
                        ph[counter2,g+1-j] = float(equilibrium_data.loc[g,'data'])
                    counter2 = counter2 + 1
                
                equilibrium_data_dict[counter] = ph
                counter = counter + 1
            ##############################################
            
            #################### IO Flow Values ########################
            counter = 1
            equilibrium_data_dict_IO = {}
            for i in range(0,(self.num_arcs*self.num_players+1)*num_od_pairs,\
                           self.num_arcs*self.num_players+1):
                equilibrium_data = pandas.read_csv(IO_cost_flow_values[iteration],\
                                                   skiprows=i,nrows=self.num_arcs*self.num_players+1)
                #pdb.set_trace()
                sets = []
                for j in range(0,self.num_arcs*self.num_players,self.num_arcs):
                    sets.append((j,j+self.num_arcs))
                
                ph = {}
                counter2 = 1
                for (j,k) in sets:
                    for g in range(j,k):
                        ph[counter2,g+1-j] = float(equilibrium_data.loc[g,'data'])
                    counter2 = counter2 + 1
                
                equilibrium_data_dict_IO[counter] = ph
                counter = counter + 1
            ##############################################
            
            ############ Calculating the Error ###################
            ph = 0
            for od in range(1,num_od_pairs+1):
                for i in range(1,self.num_players+1):
                    for j in range(1,self.num_arcs+1):
                        #pdb.set_trace()
                        ph = ph + (equilibrium_data_dict[od][i,j] - equilibrium_data_dict_IO[od][i,j])**2
            
            ph = math.sqrt(ph)
            
            print("difference",str(iteration),ph)
            
            differences[iteration] = ph
        
        print("FINAL DIFFERENCES", differences)
        
        ############### Saving the Differences #############
        self.flow_differences = differences
    
    def creating_boxplots(self,if_different_costs=0):
        # Function doesn't create boxplots any more, instead creates csv files with data
        # The if_different_costs parameter is also no longer relevant
        
        ############ Saving the Values to a CSV File ###########
        #https://www.geeksforgeeks.org/writing-csv-files-in-python/
        
        ########## Parameter dictionary ############
        
        #########################################
        collective_dict = [self.errors_c,self.errors_chat,self.flow_differences,self.timing,self.obj_func_values,self.return_values]
        headers = [i for i in range(1,self.num_trials+1)]
        
        name_of_file = "data_saving_experiment_different_costs_"+self.name_of_graph+"_num_players_"\
                    +str(self.num_players)+"_alpha_"+str(round(self.alpha))\
                    +".csv"
        
        with open(name_of_file,'w') as our_file:
            csv_house = csv.DictWriter(our_file,fieldnames=headers)
            csv_house.writeheader()
            csv_house.writerows(collective_dict)