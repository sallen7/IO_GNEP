####### Class for RK Code ###########

import pdb #for debugging
import numpy as np                     
import pyomo.environ as pyo
from pyomo.opt import SolverFactory 
from pyomo.opt import SolverStatus, TerminationCondition
import pyomo.mpec as pyompec #for the complementarity
import math


class RK_IO_methods():
    
    def __init__(self,data={},amount_of_data=100,number_of_nodes=16,\
                 number_of_arcs=48,incidence_matrix=np.zeros((16,48)),\
                 lowerbound=0,upperbound=1,lowerbound_c_hat=2,\
                 upperbound_c_hat=10,alpha=1.5):
        ## The data will be assumed to be a giant dictionary of dictionaries
        self.data = data
        self.amount_of_data = amount_of_data
        self.number_of_nodes = number_of_nodes
        self.number_of_arcs = number_of_arcs
        self.incidence_matrix = incidence_matrix
        self.lowerbound = lowerbound
        self.upperbound = upperbound
        self.lowerbound_c_hat = lowerbound_c_hat
        self.upperbound_c_hat = upperbound_c_hat
        self.alpha = alpha
    
    ###########################################################################
    ################# I. Methods that Set up the Models #######################
    ###########################################################################
    
    def RK_model_n_players(self,n=5):
        #n represents the number of players
        
        #For n players, the data has to be formatted as (,) with the 
        #row numbers as the players and the column numbers as the arc numbers
        
        #The c is still the same across all players
        
        RK_model = pyo.ConcreteModel()
        RK_model.data_index = pyo.RangeSet(1,self.amount_of_data)
        
        def block_rule(model,sample):
            model.node_index = pyo.RangeSet(1,self.number_of_nodes) #number of nodes
            model.arc_index = pyo.RangeSet(1,self.number_of_arcs) #number of arcs
            model.player_index = pyo.RangeSet(1,n) #number of players
            
            model.x = pyo.Param(model.player_index,model.arc_index,initialize=self.data[sample])
            model.v = pyo.Var(model.player_index,model.node_index)
            model.u_players = pyo.Var(model.player_index,model.arc_index,domain=pyo.NonNegativeReals)
            model.u = pyo.Var(model.arc_index,domain=pyo.NonNegativeReals)
            model.c = pyo.Var(model.arc_index)
            model.c_hat = pyo.Var(model.arc_index) #same across all n players
            
            ######### Data - Importing in the Incidence Mat ############
            def augmented_incidence_mat_func(block_model,i,j): #need both the node and var indices
                return int(self.incidence_matrix[i-1,j-1])

            model.incidence_mat = pyo.Param(model.node_index,\
                        model.arc_index,rule=augmented_incidence_mat_func)
            
            ########## Stationarity Condition ####################
            def stationary_condition(block_model,p,i): #p represents the players, i represents the arcs
                return block_model.c[i]*sum(block_model.x[p1,i] for p1 in range(1,n+1)) \
                + block_model.c[i]*block_model.x[p,i] + block_model.c_hat[i] +\
                + sum(block_model.incidence_mat[j,i]*block_model.v[p,j] for j in range(1,self.number_of_nodes+1))\
                        - block_model.u_players[p,i] + block_model.u[i]
            
            model.stationary_expression = pyo.Expression(model.player_index,model.arc_index,\
                                                         rule=stationary_condition)
            
            model.t1 = pyo.Var(model.player_index,model.arc_index,domain=pyo.NonNegativeReals)
            
            def upperbound_t1(block_model,p,i):
                return block_model.stationary_expression[p,i] <= block_model.t1[p,i]
            
            def lowerbound_t1(block_model,p,i):
                return -1*block_model.t1[p,i] <= block_model.stationary_expression[p,i]
            
            model.t1_upperbound = pyo.Constraint(model.player_index,model.arc_index,\
                                                 rule=upperbound_t1)
            
            model.t1_lowerbound = pyo.Constraint(model.player_index,model.arc_index,\
                                                 rule=lowerbound_t1)
            
            ########## Complementarity Conditions #################
            def comp_condition(block_model,p,i): 
                return block_model.x[p,i]*block_model.u_players[p,i]
            
            model.comp_expression = pyo.Expression(model.player_index,model.arc_index,\
                                                   rule=comp_condition)
            
            model.t2 = pyo.Var(model.player_index,model.arc_index,domain=pyo.NonNegativeReals)
            
            def upperbound_t2(block_model,p,i):
                return block_model.comp_expression[p,i] <= block_model.t2[p,i]
            
            def lowerbound_t2(block_model,p,i):
                return -1*block_model.t2[p,i] <= block_model.comp_expression[p,i]

            model.t2_upperbound = pyo.Constraint(model.player_index,model.arc_index,\
                                                 rule=upperbound_t2)
            
            model.t2_lowerbound = pyo.Constraint(model.player_index,model.arc_index,\
                                                 rule=lowerbound_t2)
    
            
            ############# Joint Constraint ########################
            def joint_condition(block_model,i):
                return (self.alpha - sum(block_model.x[p1,i] for p1 in range(1,n+1)))*(block_model.u[i])
            
            model.joint_constraint = pyo.Expression(model.arc_index,rule=joint_condition)
            
            model.t3 = pyo.Var(model.arc_index,domain=pyo.NonNegativeReals)
            
            def upperbound_t3(block_model,i):
                return block_model.joint_constraint[i] <= block_model.t3[i]
            
            def lowerbound_t3(block_model,i):
                return -1*block_model.t3[i] <= block_model.joint_constraint[i]
            
            model.t3_upperbound = pyo.Constraint(model.arc_index,rule=upperbound_t3)
            model.t3_lowerbound = pyo.Constraint(model.arc_index,rule=lowerbound_t3)
            
        #############################################################################
        ####################### End of Block Rule ###################################
        #############################################################################
        
        RK_model.blocks_of_residuals = pyo.Block(RK_model.data_index,rule=block_rule)
        
        ################## Linking the c variables AND Add bounds #######################
        def linking_c_variables(model,i,j):
            if i < self.amount_of_data:
                return model.blocks_of_residuals[i].c[j] == model.blocks_of_residuals[i+1].c[j]
            else:
                return pyo.Constraint.Skip
            
        RK_model.linking_c = pyo.Constraint(RK_model.data_index,\
                                             RK_model.blocks_of_residuals[1].arc_index,\
                                             rule=linking_c_variables)
        
        def upperbound_c(model,j):
            return model.blocks_of_residuals[1].c[j] <= self.upperbound
        
        def lowerbound_c(model,j):
            return self.lowerbound <= model.blocks_of_residuals[1].c[j] 
        
        RK_model.c_lower_bound = pyo.Constraint(RK_model.blocks_of_residuals[1].arc_index,\
                                                 rule=lowerbound_c)
        
        RK_model.c_upper_bound = pyo.Constraint(RK_model.blocks_of_residuals[1].arc_index,\
                                                 rule=upperbound_c)
        
        ################## Linking the c_hat variables AND Add bounds #######################
        def linking_c_hat_variables(model,i,j):
            if i < self.amount_of_data:
                return model.blocks_of_residuals[i].c_hat[j] == model.blocks_of_residuals[i+1].c_hat[j]
            else:
                return pyo.Constraint.Skip
            
        RK_model.linking_c_hat = pyo.Constraint(RK_model.data_index,\
                                             RK_model.blocks_of_residuals[1].arc_index,\
                                             rule=linking_c_hat_variables)
        
        def upperbound_c_hat(model,j):
            return model.blocks_of_residuals[1].c_hat[j] <= self.upperbound_c_hat
        
        def lowerbound_c_hat(model,j):
            return self.lowerbound_c_hat <= model.blocks_of_residuals[1].c_hat[j] 
        
        RK_model.c_lower_bound_hat = pyo.Constraint(RK_model.blocks_of_residuals[1].arc_index,\
                                                 rule=lowerbound_c_hat)
        
        RK_model.c_upper_bound_hat = pyo.Constraint(RK_model.blocks_of_residuals[1].arc_index,\
                                                 rule=upperbound_c_hat)
        
        ################## Objective Function ################################
        
        def objective_function(model):
            return sum(sum(sum(model.blocks_of_residuals[sample].t1[p1,i]\
            for sample in range(1,self.amount_of_data+1)) for i in range(1,self.number_of_arcs+1)) \
            for p1 in range(1,n+1)) +\
            sum(sum(sum(model.blocks_of_residuals[sample].t2[p1,i]\
            for sample in range(1,self.amount_of_data+1)) for i in range(1,self.number_of_arcs+1))\
            for p1 in range(1,n+1)) +\
            sum(sum(model.blocks_of_residuals[sample].t3[i]\
            for sample in range(1,self.amount_of_data+1)) for i in range(1,self.number_of_arcs+1))
            
        RK_model.obj_func_expression = pyo.Expression(rule=objective_function)
        
        RK_model.obj_func = pyo.Objective(expr=RK_model.obj_func_expression)
        
        self.RK_model = RK_model.clone()
        
    def RK_model_n_players_different_costs(self,n=5):
        #n represents the number of players
        
        #For n players, the data has to be formatted as (,) with the 
        #row numbers as the players and the column numbers as the arc numbers
        
        #The c is still the same across all players
        
        RK_model = pyo.ConcreteModel()
        RK_model.data_index = pyo.RangeSet(1,self.amount_of_data)
        
        def block_rule(model,sample):
            model.node_index = pyo.RangeSet(1,self.number_of_nodes) #number of nodes
            model.arc_index = pyo.RangeSet(1,self.number_of_arcs) #number of arcs
            model.player_index = pyo.RangeSet(1,n) #number of players
            
            model.x = pyo.Param(model.player_index,model.arc_index,initialize=self.data[sample])
            model.v = pyo.Var(model.player_index,model.node_index)
            model.u_players = pyo.Var(model.player_index,model.arc_index,domain=pyo.NonNegativeReals)
            model.u = pyo.Var(model.arc_index,domain=pyo.NonNegativeReals)
            model.c = pyo.Var(model.player_index,model.arc_index)
            model.c_hat = pyo.Var(model.player_index,model.arc_index)
            
            ######### Data - Importing in the Incidence Mat ############
            def augmented_incidence_mat_func(block_model,i,j): #need both the node and var indices
                return int(self.incidence_matrix[i-1,j-1])

            model.incidence_mat = pyo.Param(model.node_index,\
                        model.arc_index,rule=augmented_incidence_mat_func)
            
            ############################ Stationarity Condition ############################
            def stationary_condition(block_model,p,i): #p represents the players, i represents the arcs
                return block_model.c[p,i]*sum(block_model.x[p1,i] for p1 in range(1,n+1)) \
                + block_model.c[p,i]*block_model.x[p,i] \
                + sum(block_model.incidence_mat[j,i]*block_model.v[p,j] for j in range(1,self.number_of_nodes+1))\
                        + block_model.c_hat[p,i] - block_model.u_players[p,i] + block_model.u[i]
            
            model.stationary_expression = pyo.Expression(model.player_index,model.arc_index,\
                                                         rule=stationary_condition)
            
            ######## L1 Norm for the the Stationarity Condition ##########
            model.t1 = pyo.Var(model.player_index,model.arc_index,domain=pyo.NonNegativeReals)
            
            def upper_bound_t1(block_model,p,i):
                return block_model.stationary_expression[p,i] <= block_model.t1[p,i]
            
            def lower_bound_t1(block_model,p,i):
                return -1*block_model.t1[p,i] <= block_model.stationary_expression[p,i] 
            
            model.upper_bound_t1_constraint = pyo.Constraint(model.player_index,\
                                                             model.arc_index,\
                                                             rule=upper_bound_t1)
            
            model.lower_bound_t1_constraint = pyo.Constraint(model.player_index,\
                                                             model.arc_index,\
                                                             rule=lower_bound_t1)           
            
            ########## Complementarity Conditions #################
            def comp_condition(block_model,p,i):
                return block_model.x[p,i]*block_model.u_players[p,i]
            
            model.comp_expression = pyo.Expression(model.player_index,model.arc_index,\
                                                   rule=comp_condition)
            
            model.t2 = pyo.Var(model.player_index,model.arc_index,domain=pyo.NonNegativeReals)
            
            def upper_bound_t2(block_model,p,i):
                return block_model.comp_expression[p,i] <= block_model.t2[p,i]
            
            def lower_bound_t2(block_model,p,i):
                return -1*block_model.t2[p,i] <= block_model.comp_expression[p,i] 
            
            model.upper_bound_t2_constraint = pyo.Constraint(model.player_index,\
                                                             model.arc_index,\
                                                             rule=upper_bound_t2)
            
            model.lower_bound_t2_constraint = pyo.Constraint(model.player_index,\
                                                             model.arc_index,\
                                                             rule=lower_bound_t2)
            
            def joint_condition(block_model,i):
                return (self.alpha - sum(block_model.x[p1,i] for p1 in range(1,n+1)))*(block_model.u[i])
            
            model.joint_constraint = pyo.Expression(model.arc_index,rule=joint_condition)
            
            model.t3 = pyo.Var(model.arc_index,domain=pyo.NonNegativeReals)
            
            def upper_bound_t3(block_model,i):
                return block_model.joint_constraint[i] <= block_model.t3[i]
            
            def lower_bound_t3(block_model,i):
                return -1*block_model.t3[i] <= block_model.joint_constraint[i]
            
            model.upper_bound_t3_constraint = pyo.Constraint(model.arc_index,\
                                                             rule=upper_bound_t3)
            
            model.lower_bound_t3_constraint = pyo.Constraint(model.arc_index,\
                                                             rule=lower_bound_t3)
            
        #############################################################################
        ####################### End of Block Rule ###################################
        #############################################################################
        
        RK_model.blocks_of_residuals = pyo.Block(RK_model.data_index,rule=block_rule)
        
        ################## Linking the c variables AND Add bounds #######################
        def linking_c_variables(model,p,i,j):
            if i < self.amount_of_data:
                return model.blocks_of_residuals[i].c[p,j] == model.blocks_of_residuals[i+1].c[p,j]
            else:
                return pyo.Constraint.Skip
            
        RK_model.linking_c = pyo.Constraint(RK_model.blocks_of_residuals[1].player_index,\
                                            RK_model.data_index,\
                                             RK_model.blocks_of_residuals[1].arc_index,\
                                             rule=linking_c_variables)
        
        def upperbound_c(model,p,j):
            return model.blocks_of_residuals[1].c[p,j] <= self.upperbound
        
        def lowerbound_c(model,p,j):
            return self.lowerbound <= model.blocks_of_residuals[1].c[p,j] 
        
        RK_model.c_lower_bound = pyo.Constraint(RK_model.blocks_of_residuals[1].player_index,\
                                                RK_model.blocks_of_residuals[1].arc_index,\
                                                 rule=lowerbound_c)
        
        RK_model.c_upper_bound = pyo.Constraint(RK_model.blocks_of_residuals[1].player_index,\
                                                RK_model.blocks_of_residuals[1].arc_index,\
                                                 rule=upperbound_c)
        
        ################### Linking c_hat and Adding Bounds ##################
        def linking_c_hat_variables(model,p,i,j):
            if i < self.amount_of_data:
                return model.blocks_of_residuals[i].c_hat[p,j] == model.blocks_of_residuals[i+1].c_hat[p,j]
            else:
                return pyo.Constraint.Skip
            
        RK_model.linking_c_hat = pyo.Constraint(RK_model.blocks_of_residuals[1].player_index,\
                                            RK_model.data_index,\
                                             RK_model.blocks_of_residuals[1].arc_index,\
                                             rule=linking_c_hat_variables)
        
        def upperbound_c_hat(model,p,j):
            return model.blocks_of_residuals[1].c_hat[p,j] <= self.upperbound_c_hat
        
        def lowerbound_c_hat(model,p,j):
            return self.lowerbound_c_hat <= model.blocks_of_residuals[1].c_hat[p,j] 
        
        RK_model.c_hat_lower_bound = pyo.Constraint(RK_model.blocks_of_residuals[1].player_index,\
                                                RK_model.blocks_of_residuals[1].arc_index,\
                                                 rule=lowerbound_c_hat)
        
        RK_model.c_hat_upper_bound = pyo.Constraint(RK_model.blocks_of_residuals[1].player_index,\
                                                RK_model.blocks_of_residuals[1].arc_index,\
                                                 rule=upperbound_c_hat)
        
        ################## Objective Function ################################
        
        def objective_function(model):
            return sum(sum(sum(model.blocks_of_residuals[sample].t1[p1,i]\
            for sample in range(1,self.amount_of_data+1)) for p1 in range(1,n+1)) \
            for i in range(1,self.number_of_arcs+1)) +\
            sum(sum(sum(model.blocks_of_residuals[sample].t2[p1,i]\
            for sample in range(1,self.amount_of_data+1)) for p1 in range(1,n+1))\
            for i in range(1,self.number_of_arcs+1)) +\
            sum(sum(model.blocks_of_residuals[sample].t3[i] \
            for sample in range(1,self.amount_of_data+1)) for i in range(1,self.number_of_arcs+1))
            
        RK_model.obj_func_expression = pyo.Expression(rule=objective_function)
        
        RK_model.obj_func = pyo.Objective(expr=RK_model.obj_func_expression)
        
        self.RK_model = RK_model.clone()
    
    ###########################################################################
    ################# II. Methods that SOLVE the Models ########################
    ###########################################################################
        
    def RK_model_n_player_solve(self,epsilon=1e-12):
        RK_model = self.RK_model.clone()
        print("***********Solving****************")
        solver = SolverFactory("gurobi")
        #solver.options["tol"] = epsilon #for IPOPT solver
        solver.options["BarQCPConvTol"] = 1e-14
        solver.options["BarConvTol"] = 1e-14
        output = solver.solve(RK_model,tee=True)
        
        if output.solver.termination_condition != TerminationCondition.optimal:
            print("PROBLEM: Non optimal solution")
            pdb.set_trace() #http://www.pyomo.org/blog/2015/1/8/accessing-solver
    
        print("c values",RK_model.blocks_of_residuals[1].c.extract_values())
        print("c_hat values",RK_model.blocks_of_residuals[1].c_hat.extract_values())
        print("obj func value:",pyo.value(RK_model.obj_func_expression))
        
        self.c_values = RK_model.blocks_of_residuals[1].c.extract_values()
        self.c_hat_values = RK_model.blocks_of_residuals[1].c_hat.extract_values()
        self.obj_func_value = pyo.value(RK_model.obj_func_expression)
        self.return_value = output.solver.termination_condition
        self.RK_model_solved = RK_model.clone()
        
    def RK_model_n_player_different_costs_solve(self,epsilon=1e-12):
        RK_model = self.RK_model.clone()
        print("***********Solving****************")
        solver = SolverFactory("gurobi")
        #solver.options["tol"] = epsilon #for IPOPT solver
        solver.options["BarQCPConvTol"] = 1e-14
        solver.options["BarConvTol"] = 1e-14
        output = solver.solve(RK_model,tee=True)
        
        if output.solver.termination_condition != TerminationCondition.optimal:
            print("PROBLEM: Non optimal solution")
            pdb.set_trace() #http://www.pyomo.org/blog/2015/1/8/accessing-solver
    
        print("c values",RK_model.blocks_of_residuals[1].c.extract_values())
        print("c hat values",RK_model.blocks_of_residuals[1].c_hat.extract_values())
        
        self.c_matrix = RK_model.blocks_of_residuals[1].c.extract_values()
        self.c_hat_matrix = RK_model.blocks_of_residuals[1].c_hat.extract_values()
        self.obj_func_value = pyo.value(RK_model.obj_func_expression)
        self.return_value = output.solver.termination_condition
        self.RK_model_solved = RK_model.clone()