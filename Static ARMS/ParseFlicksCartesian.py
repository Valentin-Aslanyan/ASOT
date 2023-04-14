

import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from mayavi import mlab
from scipy.interpolate import griddata


#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return data_in[:,6]#np.sqrt(data_in[:,6]**2+data_in[:,7]**2+data_in[:,8]**2)

target_z=0.0
size_x_grid=500
size_y_grid=500


time,ntblks,nlblks,newgrd,coord_x,coord_y,coord_z,data=read_flicks_file("./","flicks.002200")
new_x,new_y,new_data=z_slice(target_z,coord_x,coord_y,coord_z,data)
target_data=data_slice(new_data)


min_new_x=min(new_x)
max_new_x=max(new_x)
min_new_y=min(new_y)
max_new_y=max(new_y)


grid_x=np.linspace(min_new_x,max_new_x,num=size_x_grid)
grid_y=np.linspace(min_new_y,max_new_y,num=size_y_grid)
grid_x,grid_y=np.meshgrid(grid_x,grid_y)
grid_z=np.zeros((np.shape(grid_x)))

data_grid=griddata((new_x,new_y), target_data, (grid_x,grid_y), method='linear')


sphere_mesh=mlab.mesh(grid_x,grid_y,grid_z,scalars=data_grid[:,:],colormap='hsv')
mlab.view(azimuth=0, elevation=110, roll=90, distance=4.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)
mlab.show()


plt.pcolormesh(grid_x,grid_y,data_grid)
plt.show()





