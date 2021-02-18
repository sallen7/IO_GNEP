########## Script 3: Calculating the Total Flow Error ###################

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

############ Step 1: Importing Object ##########################
#https://www.datacamp.com/community/tutorials/pickle-python-tutorial
file_to_be_read = open("class_object_2","rb")
generalized_framework_object = pickle.load(file_to_be_read)
file_to_be_read.close()

######### Step 2: Calculate the Total Flow Error and Create Graphs ##############

generalized_framework_object.calculating_flow_error()
generalized_framework_object.creating_boxplots(if_different_costs=generalized_framework_object.if_different_costs)
#^function doesn't create boxplots any more, instead creates csv files with data

########### Step 3: Saving the Object Again ###################
#https://www.datacamp.com/community/tutorials/pickle-python-tutorial
name_of_file = "class_object_3"
test = open(name_of_file,'wb')
pickle.dump(generalized_framework_object,test)
test.close()

