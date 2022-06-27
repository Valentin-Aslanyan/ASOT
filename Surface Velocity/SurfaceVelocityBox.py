
"""
Display box surrounding surface flow (directly from ARMS control file) relative to Open/Closed boundary (or other "cage", e.g. Polarity Inversion Line)
"""

CNT_file="./arms.cnt"
cage_file="./SurfaceLineData/SLine.bin"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


Time_profiles,parameters=parse_sflow_CNT(CNT_file)

static_cages,dynamic_cages=surfaceline_read_binary_cage(cage_file)
for idx in range(len(static_cages)):
	plt.plot(static_cages[idx][0,:],static_cages[idx][1,:]-90.0,color="black")
for idx in range(len(dynamic_cages)):
	plt.plot(dynamic_cages[idx][0,:],dynamic_cages[idx][1,:]-90.0,color="black")

for idx in range(len(parameters[:,0])):
	if Time_profiles[idx]>0:
		theta_l=(1.0-parameters[idx,5])*180.0-90.0
		theta_r=(1.0-parameters[idx,6])*180.0-90.0
		phi_l=parameters[idx,10]*180.0
		phi_r=parameters[idx,11]*180.0
		plt.plot([phi_l,phi_r,phi_r,phi_l,phi_l],[theta_l,theta_l,theta_r,theta_r,theta_l])


plt.show()
