import pickle
import pdb

file_to_be_read = open("class_object_3","rb")
generalized_framework_object = pickle.load(file_to_be_read)
file_to_be_read.close()

pdb.set_trace()