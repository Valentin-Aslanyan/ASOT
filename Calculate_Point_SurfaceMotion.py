


solar_Radius=1.0*7E10
ignore_uniform_velocity=True
CNT_file="./arms.cnt"
t_start=5300.0
t_end=5299.0
delta_t=-1.0	#Set negative to advect point backwards

#Flicks coordinates
theta=+1.970302
phi=-0.1319


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *


Time_profiles,sflow_parameters=parse_sflow_CNT(CNT_file)

theta_out=1.0-theta/np.pi
phi_out=phi/np.pi

t=t_start
num_t=int((t_end-t_start)/delta_t)
if num_t<=0.0:
	print("Warning, start/end/delta wrong!")
for idx_t in range(num_t):
	new_pos=surface_move_rk_4(solar_Radius,theta_out,phi_out,t,delta_t,Time_profiles,sflow_parameters)
	theta_out=new_pos[0]
	phi_out=new_pos[1]
	t+=delta_t
	print(t)

if abs(t_end-t)>0.0:
	new_pos=surface_move_rk_4(solar_Radius,theta_out,phi_out,t,t_end-t,Time_profiles,sflow_parameters)
	theta_out=new_pos[0]
	phi_out=new_pos[1]

print(np.pi*phi_out)
print(np.pi*(1.0-theta_out))
	





