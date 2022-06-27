
"""
Advect dynamic "cages" according to ideal surface flow velocity (defined in ARMS control file)
"""

CNT_file="./arms.cnt"
SurfaceLine_initial="./SurfaceLineData/SLine.bin"
t_start=0.0
delta_t=10.0
solar_Radius=1.0*7E10
plot_t=[1000.0,2000.0,3000.0,4000.0,5000.0,6000.0,7000.0,8000.0] #Make sure times are in order and correspond correctly to filenames
SurfaceLine_outputs=[
"./SurfaceLineData/SLine_1000.bin",
"./SurfaceLineData/SLine_2000.bin",
"./SurfaceLineData/SLine_3000.bin",
"./SurfaceLineData/SLine_4000.bin",
"./SurfaceLineData/SLine_5000.bin",
"./SurfaceLineData/SLine_6000.bin",
"./SurfaceLineData/SLine_7000.bin",
"./SurfaceLineData/SLine_8000.bin"]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


static_cages,dynamic_cages_start=surfaceline_read_binary_cage(SurfaceLine_initial)
dynamic_cages=[]
for idx in range(len(dynamic_cages_start)):
	dynamic_cages.append(surfaceline_advance_cage(dynamic_cages_start[idx],t_start,plot_t,delta_t,CNT_file,solar_Radius))

for idx in range(min(len(plot_t),len(SurfaceLine_outputs))):
	dynamic_cages_temp=[]
	for idx2 in range(len(dynamic_cages_start)):
		dynamic_cages_temp.append(dynamic_cages[idx2][:,:,idx])
	surfaceline_write_binary_cage(SurfaceLine_outputs[idx],static_cages,dynamic_cages_temp)






