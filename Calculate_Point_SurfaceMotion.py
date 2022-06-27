
"""
Track positions of points on the photosphere, advected by applied surface velocity
Use this to track points which are known to be e.g. interchange reconnected at the end
"""

solar_Radius=1.0*7E10
ignore_uniform_velocity=True
CNT_file="./arms.cnt"
t_start=1236000.0
t_end=1200000.0
delta_t=-10.0	#Set negative to advect point backwards

coordinates_input='carrington'	#'flicks' - radians; 'carrington' - degrees, theta in -90,90; 'control' - same format as flows in arms.cnt
theta=[+26.828,+26.425,+26.165,+26.157,+26.591,  +24.96,+24.62,+24.08,+23.43,+21.91]
phi=[-9.7230,-9.2380,-8.6530,-8.1200,-7.5430,  +2.572,+1.671,+0.632,-0.126,-0.351]
coordinates_output='carrington'


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *


Time_profiles,sflow_parameters=parse_sflow_CNT(CNT_file)

theta_out,phi_out=change_angular_coords(np.array(theta),np.array(phi),coordinates_input,'control')

t=t_start
num_t=int((t_end-t_start)/delta_t)
if num_t<=0.0:
	print("Warning, start/end/delta wrong!")
for idx_t in range(num_t):
	for idx in range(len(theta_out)):
		new_pos=surface_move_rk_4(solar_Radius,theta_out[idx],phi_out[idx],t,delta_t,Time_profiles,sflow_parameters)
		theta_out[idx]=new_pos[0]
		phi_out[idx]=new_pos[1]
	t+=delta_t
	#print(t)

if abs(t_end-t)>0.0:
	for idx in range(len(theta_out)):
		new_pos=surface_move_rk_4(solar_Radius,theta_out[idx],phi_out[idx],t,t_end-t,Time_profiles,sflow_parameters)
		theta_out[idx]=new_pos[0]
		phi_out[idx]=new_pos[1]

theta_out,phi_out=change_angular_coords(theta_out,phi_out,'control',coordinates_output)
print("theta =",theta_out)
print("phi =",phi_out)
	





