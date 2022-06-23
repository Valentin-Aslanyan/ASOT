
"""
ARMS quantity at constant azimuthal angle (phi) on a plane (sector of annulus), optionally with fieldlines
"""

import numpy as np
target_phi=0.0	#degrees
file_directory="./"
flicks_file="flicks.0000000"
theta_lineout=10.0	#None or number in degrees, 
R_limits=None		#None or Solar radii

#Carrington coordinates
plot_gridfieldlines=True
gridfieldlines_lims=[1.0,3.0]
R_start    =[ +1.000, +1.000, +1.000, +1.000]
theta_start=[+11.471]
phi_start  =[ +0.974]

#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return np.sqrt(data_in[:,6]**2+data_in[:,7]**2+data_in[:,8]**2)


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from scipy.interpolate import griddata


size_x_grid=500
size_y_grid=500

time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)

new_logR,new_theta,new_data=phi_slice(target_phi,coord_logR,coord_theta,coord_phi,data)

if R_limits!=None:
	new_logR_full=np.copy(new_logR)
	new_theta_full=np.copy(new_theta)
	new_data_full=np.copy(new_data)
	idx_limited=0
	for idx in range(len(new_logR_full)):
		if np.exp(new_logR_full[idx])>=R_limits[0] and np.exp(new_logR_full[idx])<=R_limits[1]:
			new_logR[idx_limited]=new_logR_full[idx]
			new_theta[idx_limited]=new_theta_full[idx]
			new_data[idx_limited,:]=new_data_full[idx,:]
			idx_limited+=1
	new_logR=new_logR[:idx_limited]
	new_theta=new_theta[:idx_limited]
	new_data=new_data[:idx_limited,:]


target_x=np.exp(new_logR)*np.sin(new_theta)
target_y=np.exp(new_logR)*np.cos(np.pi-new_theta)
target_data=data_slice(new_data)

min_new_logR=min(new_logR)
max_new_logR=max(new_logR)
min_new_theta=min(new_theta)
max_new_theta=max(new_theta)
min_newR2=np.exp(min_new_logR)**2
max_new_R=np.exp(max_new_logR)


grid_x=np.linspace(min(target_x),max(target_x),num=size_x_grid)
grid_y=np.linspace(min(target_y),max(target_y),num=size_y_grid)
grid_x,grid_y=np.meshgrid(grid_x,grid_y)
data_grid=griddata((target_x,target_y), target_data, (grid_x,grid_y), method='linear')
for idx_x in range(len(data_grid[0,:])):
	for idx_y in range(len(data_grid[:,0])):
		if grid_x[0,idx_x]**2+grid_y[idx_y,0]**2<min_newR2:
			data_grid[idx_y,idx_x]=np.nan

num_arc=200
arc_r=np.concatenate((np.exp(min_new_logR)*np.ones((num_arc)),np.exp(max_new_logR)*np.ones((num_arc)),[np.exp(min_new_logR)]))
arc_theta=np.concatenate((np.linspace(min_new_theta,max_new_theta,num=num_arc),np.linspace(max_new_theta,min_new_theta,num=num_arc),[min_new_theta]))
arc_x=arc_r*np.sin(arc_theta)
arc_y=arc_r*np.cos(np.pi-arc_theta)


if plot_gridfieldlines:
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	B_flicks=data[:,:,:,:,6:9]
	B_flicks[:,:,:,:,1]*=-1.0
	plt.figure("ARMS grid and field lines",figsize=(3.1*1.5,6*1.5))
	for idx in range(len(coord_logR)):
		if is_block_in_limits(coord_logR[idx],coord_theta[idx],coord_phi[idx],target_phi,R_limits):
			x1=np.exp(coord_logR[idx][0])*np.sin(coord_theta[idx][0])
			z1=np.exp(coord_logR[idx][0])*np.cos(np.pi-coord_theta[idx][0])
			x2=np.exp(coord_logR[idx][1])*np.sin(coord_theta[idx][0])
			z2=np.exp(coord_logR[idx][1])*np.cos(np.pi-coord_theta[idx][0])
			x3=np.exp(coord_logR[idx][1])*np.sin(coord_theta[idx][1])
			z3=np.exp(coord_logR[idx][1])*np.cos(np.pi-coord_theta[idx][1])
			x4=np.exp(coord_logR[idx][0])*np.sin(coord_theta[idx][1])
			z4=np.exp(coord_logR[idx][0])*np.cos(np.pi-coord_theta[idx][1])
			plt.plot([x1,x2],[z1,z2],color="grey",linewidth=1)
			plt.plot([x3,x4],[z3,z4],color="grey",linewidth=1)
			anglearc=np.linspace(coord_theta[idx][0],coord_theta[idx][1],num=100)
			arc1x=np.exp(coord_logR[idx][0])*np.sin(anglearc)
			arc1y=np.exp(coord_logR[idx][0])*np.cos(np.pi-anglearc)
			arc2x=np.exp(coord_logR[idx][1])*np.sin(anglearc)
			arc2y=np.exp(coord_logR[idx][1])*np.cos(np.pi-anglearc)
			plt.plot(arc1x,arc1y,color="grey",linewidth=1)
			plt.plot(arc2x,arc2y,color="grey",linewidth=1)
	for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
		th_temp,ph_temp=change_angular_coords(theta_start[idx],phi_start[idx],from_type='carrington',to_type='flicks')
		field_line_start=np.array([R_start[idx],th_temp+0.5*np.pi,ph_temp])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
		plt.plot(field_line_X,field_line_Z,color="red",linewidth=2)
	plt.axis('off')
	plt.xlim([0,max_new_R*1.05])
	plt.ylim([-max_new_R,max_new_R])
	plt.savefig("Grid_Field_XZ.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


if theta_lineout!=None:
	R_lineout=np.linspace(np.exp(min(new_logR)),np.exp(max(new_logR)),num=100)
	data_lineout=griddata((new_logR,new_theta), target_data, (np.log(R_lineout),(90.0-theta_lineout)*DEG2RAD*np.ones((len(R_lineout)))), method='linear')
	plt.figure("Theta lineout",figsize=(6.0,6.0))
	plt.plot(R_lineout,data_lineout)
	plt.savefig("ThetaLineoutPhi2D.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


plt.figure("data plot",figsize=(3.1*1.5,6*1.5))
plt.pcolormesh(grid_x,grid_y,data_grid,cmap='RdPu',rasterized=True)#,vmin=0.0,vmax=1.0)
if plot_gridfieldlines:
	for idx in range(len(coord_logR)):
		#if coord_phi[idx][0]<=target_phi/180.0*np.pi and coord_phi[idx][1]>=target_phi/180.0*np.pi:
		if is_block_in_limits(coord_logR[idx],coord_theta[idx],coord_phi[idx],target_phi,R_limits):
			x1=np.exp(coord_logR[idx][0])*np.sin(coord_theta[idx][0])
			z1=np.exp(coord_logR[idx][0])*np.cos(np.pi-coord_theta[idx][0])
			x2=np.exp(coord_logR[idx][1])*np.sin(coord_theta[idx][0])
			z2=np.exp(coord_logR[idx][1])*np.cos(np.pi-coord_theta[idx][0])
			x3=np.exp(coord_logR[idx][1])*np.sin(coord_theta[idx][1])
			z3=np.exp(coord_logR[idx][1])*np.cos(np.pi-coord_theta[idx][1])
			x4=np.exp(coord_logR[idx][0])*np.sin(coord_theta[idx][1])
			z4=np.exp(coord_logR[idx][0])*np.cos(np.pi-coord_theta[idx][1])
			plt.plot([x1,x2],[z1,z2],color="grey",linewidth=1)
			plt.plot([x3,x4],[z3,z4],color="grey",linewidth=1)
			anglearc=np.linspace(coord_theta[idx][0],coord_theta[idx][1],num=100)
			arc1x=np.exp(coord_logR[idx][0])*np.sin(anglearc)
			arc1y=np.exp(coord_logR[idx][0])*np.cos(np.pi-anglearc)
			arc2x=np.exp(coord_logR[idx][1])*np.sin(anglearc)
			arc2y=np.exp(coord_logR[idx][1])*np.cos(np.pi-anglearc)
			plt.plot(arc1x,arc1y,color="grey",linewidth=1)
			plt.plot(arc2x,arc2y,color="grey",linewidth=1)
	for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
		th_temp,ph_temp=change_angular_coords(theta_start[idx],phi_start[idx],from_type='carrington',to_type='flicks')
		field_line_start=np.array([R_start[idx],th_temp+0.5*np.pi,ph_temp])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
		plt.plot(field_line_X,field_line_Z,color="red",linewidth=2)
plt.plot(arc_x,arc_y,color="blue",linewidth=2)
plt.axis('off')
plt.xlim([0,3.1])
plt.ylim([-3,3])
plt.savefig("FlicksPhi2D.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


plt.show()



