
"""
ARMS quantity at constant radius mapped onto a 3D sphere, with optional fieldlines
"""

import numpy as np
target_R=1.0
file_directory="./"
flicks_file="flicks.0000000"

#Carrington coordinates
plot_field_lines=True  
R_start    =[ +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000, +1.000]
theta_start=[+53.812,+26.196,+28.258,+24.019,+22.357,+25.279,+20.409,+16.627,+34.962,+62.349,+64.068,+56.849,+55.073,+54.042,+49.057,+40.004,+24.706]
phi_start  =[-46.467,+65.031,-16.157, -0.859,+12.376, -8.881,+24.351,+44.175,-29.622,-37.815, -8.766,+14.610,+22.976,+39.935,+55.462,+65.546, -9.282]


#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return np.sqrt(data_in[:,6]**2+data_in[:,7]**2+data_in[:,8]**2)


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from mayavi import mlab
from scipy.interpolate import griddata


size_t_grid=400
size_p_grid=500

time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)
if plot_field_lines:
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	B_flicks=data[:,:,:,:,6:9]
	B_flicks[:,:,:,:,1]*=-1.0

new_theta,new_phi,new_data=R_slice(target_R,coord_logR,coord_theta,coord_phi,data)
target_data=data_slice(new_data)


min_new_theta=min(new_theta)
max_new_theta=max(new_theta)
min_new_phi=min(new_phi)
max_new_phi=max(new_phi)


grid_theta=np.linspace(min_new_theta,max_new_theta,num=size_t_grid)
grid_phi=np.linspace(min_new_phi,max_new_phi,num=size_p_grid)
grid_theta,grid_phi=np.meshgrid(grid_theta,grid_phi)

data_grid=griddata((new_theta,new_phi), target_data, (grid_theta,grid_phi), method='linear')

X_grid=target_R*np.sin(grid_theta)*np.cos(grid_phi)
Y_grid=target_R*np.sin(grid_theta)*np.sin(grid_phi)
Z_grid=target_R*np.cos(grid_theta)

mlab.figure("|B| field")
sphere_mesh=mlab.mesh(X_grid[::-1,:],Y_grid[::-1,:],Z_grid,scalars=np.log(data_grid[:,:]),colormap='bone')
sphere_mesh.actor.property.backface_culling = True
sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True

if plot_field_lines:
	for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
		th_temp,ph_temp=change_angular_coords(theta_start[idx],phi_start[idx],from_type='carrington',to_type='flicks')
		field_line_start=np.array([R_start[idx],th_temp+0.5*np.pi,ph_temp])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,1.0,2.9,nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2])
		field_line_Y=-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2])
		field_line_Z=field_line_sph[:,0]*np.cos(field_line_sph[:,1])
		mlab.plot3d(field_line_X,field_line_Y,field_line_Z,line_width=0.01,color=(0,0,0),tube_radius=0.004)
		plt.plot(field_line_sph[:,2],field_line_sph[:,1],color="black")

mlab.view(azimuth=0, elevation=110, roll=90, distance=4.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)
mlab.show()



