
"""
Read in a flicks file and write a new one with a new, user-created variable and optionally change some existing variables
"""

input_directory="Start"
output_directory="End"
flicks_labels=["0000000"]
new_variable_name="1JoverB"	#Lead with integer denoting number of entries, e.g. "1Alfven Speed"

#Function to generate new variable from old arrays; see below for how to do this more generally (e.g. using adjacent blocks)
def new_variable_function(data_in):
	return np.sqrt((data_in[...,9]**2+data_in[...,10]**2+data_in[...,11]**2)/(data_in[...,6]**2+data_in[...,7]**2+data_in[...,8]**2))


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *


duplicate_flicks_header(input_directory,output_directory,new_variable_name)
for label in flicks_labels:
	time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data_old=read_flicks_file(input_directory,"flicks."+str(label))
	nvar_old=np.shape(data_old)[-1]
	data_new=np.concatenate((data_old,np.expand_dims(new_variable_function(data_old),axis=4)),axis=4)
	duplicate_flicks_file(os.path.join(input_directory,"flicks."+str(label)),os.path.join(output_directory,"flicks."+str(label)),nvar_old,data_new)




