

file_directory="./"
flicks_file="flicks.0000000"
target_z=0.0 #Take slice at this z, ARMS coordinates

import numpy as np
#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return np.sqrt(data_in[:,5]**2+data_in[:,6]**2+data_in[:,7]**2)

plot_field_lines=True  
x_start=[-2.500,-1.600,-3.700,-2.870,-3.100,-2.230,+2.500,+1.600,+3.700,+2.870,+3.100,+2.230]
y_start=[+4.000,+3.360,+1.650,+1.420,-1.970,-1.370,+4.000,+3.360,+1.650,+1.420,-1.970,-1.370]
z_start=[+0.000,+0.000,+0.000,+0.000,+0.000,+0.000,+0.000,+0.000,+0.000,+0.000,+0.000,+0.000]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from mayavi import mlab
from scipy.interpolate import griddata


size_x_grid=500
size_y_grid=500


time,ntblks,nlblks,newgrd,coord_x,coord_y,coord_z,data=read_flicks_file(file_directory,flicks_file)
print(min(coord_x.flatten()),max(coord_x.flatten()))
print(min(coord_y.flatten()),max(coord_y.flatten()))
print(min(coord_z.flatten()),max(coord_z.flatten()))
if plot_field_lines:
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	B_flicks=data[:,:,:,:,5:8]
new_x,new_y,new_data=z_slice(target_z,coord_x,coord_y,coord_z,data)
target_data=data_slice(new_data)


min_new_x=min(new_x)
max_new_x=max(new_x)
min_new_y=min(new_y)
max_new_y=max(new_y)


grid_x=np.linspace(min_new_x,max_new_x,num=size_x_grid)
grid_y=np.linspace(min_new_y,max_new_y,num=size_y_grid)
grid_x,grid_y=np.meshgrid(grid_x,grid_y)
grid_z=target_z*np.ones((np.shape(grid_x)))

data_grid=griddata((new_x,new_y), target_data, (grid_x,grid_y), method='linear')

mlab.figure(bgcolor=(1,1,1),size=(1200,1000))

sphere_mesh=mlab.mesh(grid_x,grid_y,grid_z,scalars=data_grid[:,:],colormap='hsv')
plt.pcolormesh(grid_x,grid_y,data_grid,cmap='Reds',rasterized=True)

if plot_field_lines:
	for idx in range(min(len(x_start),len(y_start),len(z_start))):
		field_line_start=np.array([x_start[idx],y_start[idx],z_start[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_x,coord_y,coord_z,-B_flicks,-1.0,100.0,nlblks,n1pm1,n2pm1,n3pm1,max_steps=2E4,step_size=1E-1,coord_type="Cart")
		field_line_X=field_line_sph[:,0]
		field_line_Y=field_line_sph[:,1]
		field_line_Z=field_line_sph[:,2]
		mlab.plot3d(field_line_X,field_line_Y,field_line_Z,line_width=0.01,color=(0,0,0),tube_radius=0.01)
		plt.plot(field_line_sph[:,0],field_line_sph[:,1],color="black")

mlab.view(azimuth=0, elevation=30, roll=0, distance=24.0, focalpoint=None)#, roll=None, reset_roll=True, figure=None)

plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$y$",fontsize=20)
plt.xlabel(r"$x$",fontsize=20)
plt.savefig("FlicksCartesian.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

mlab.show()
plt.show()





