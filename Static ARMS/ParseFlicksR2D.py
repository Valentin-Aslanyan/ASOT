
"""
ARMS quantity at constant radius displayed as a 2D image
"""

import numpy as np
target_R=1.0
file_directory="./"
flicks_file="flicks.0000000"

#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return np.sqrt(data_in[:,6]**2+data_in[:,7]**2+data_in[:,8]**2)

plot_nulls=False
nulls_directory="./"


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
if plot_nulls:
	import HQVseg
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from scipy.interpolate import griddata


size_t_grid=1600
size_p_grid=800

time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)
if plot_nulls:
	null_locations=HQVseg.nulls_from_dir(nulls_directory)

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

plt.figure("ARMS grid",figsize=(20,10))
for idx in range(len(coord_logR)):
	if np.exp(coord_logR[idx][0])<=target_R and np.exp(coord_logR[idx][1])>=target_R:
		th1=coord_theta[idx][0]*RAD2DEG
		ph1=coord_phi[idx][0]*RAD2DEG
		th2=coord_theta[idx][1]*RAD2DEG
		ph2=coord_phi[idx][1]*RAD2DEG
		plt.plot([ph1,ph1,ph2,ph2,ph1],[th1,th2,th2,th1,th1],color="grey",linewidth=1)
plt.axis('off')
plt.savefig("Grid_PhiTheta.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


fig1=plt.figure("data",figsize=(10,8))
ax1=fig1.gca()
plt.title("$R="+"{:.1f}".format(target_R)+"$",fontsize=20)
color_plot1=plt.pcolormesh(grid_phi*RAD2DEG,grid_theta*RAD2DEG-90.0,data_grid[:-1,:-1],cmap='bone_r',rasterized=True)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
if plot_nulls:
	plt.plot(null_locations[:,2],null_locations[:,1],'+',color="black")
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar1=fig1.colorbar(color_plot1)
cbar1.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
#cbar1.set_label(label="$[\mathrm{a.u.}]$",fontsize=20)
plt.savefig("FlicksR2D.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


plt.show()



