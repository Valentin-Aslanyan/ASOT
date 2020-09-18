

CNT_file="./arms.cnt"
SurfaceLine_initial="./SurfaceLineData/SLine_4000.bin"
t_start=4000.0
delta_t=10.0
solar_Radius=1.0*7E10
plot_t=[4600.0,4800.0,5000.0,5200.0,5400.0,5600.0,5800.0,6000.0] #Make sure times are in order and correspond correctly to filenames
SurfaceLine_outputs=[
"./SurfaceLineData/SLine_4600.bin",
"./SurfaceLineData/SLine_4800.bin",
"./SurfaceLineData/SLine_5000.bin",
"./SurfaceLineData/SLine_5200.bin",
"./SurfaceLineData/SLine_5400.bin",
"./SurfaceLineData/SLine_5600.bin",
"./SurfaceLineData/SLine_5800.bin",
"./SurfaceLineData/SLine_6000.bin"]


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
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






