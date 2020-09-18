

import numpy as np
target_R=np.linspace(1.01,1.6,num=100)
file_directory="./"
flicks_file="flicks.0441995"

frames_per_step=10
frames_per_sec=20
pad_start_frames=10
pad_end_frames=10

plot_field_lines=True  
R_start    =[+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000]
theta_start=[+2.510,+2.028,+2.064,+1.990,+1.961,+2.012,+1.927,+1.861,+2.181,+2.659,+2.689,+2.563,+2.532,+2.514,+2.427,+2.269,+2.002]
phi_start  =[-0.811,+1.135,-0.282,-0.015,+0.216,-0.155,+0.425,+0.771,-0.517,-0.660,-0.153,+0.255,+0.401,+0.697,+0.968,+1.144,-0.162]


#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return np.sqrt((data_in[:,9]**2+data_in[:,10]**2+data_in[:,11]**2)/(data_in[:,6]**2+data_in[:,7]**2+data_in[:,8]**2))

def R_label_function(R):
	return "$R={:.2f}$".format(R)
	


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
from scipy.interpolate import griddata
from mayavi import mlab
from subprocess import call


grid_theta=np.linspace(0.0,np.pi,num=400)
grid_phi=np.linspace(-np.pi,np.pi,num=500)
grid_theta,grid_phi=np.meshgrid(grid_theta,grid_phi)
data_list=[]
min_data=None
max_data=None
time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)
if plot_field_lines:
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	B_flicks=data[:,:,:,:,6:9]
	B_flicks[:,:,:,:,1]*=-1.0
	field_line_R=[]
	field_line_X=[]
	field_line_Y=[]
	field_line_Z=[]
	for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
		field_line_start=np.array([R_start[idx],theta_start[idx],phi_start[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,1.0,2.9,nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_R.append(field_line_sph[:,0])
		field_line_X.append(field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2]))
		field_line_Y.append(-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2]))
		field_line_Z.append(field_line_sph[:,0]*np.cos(field_line_sph[:,1]))
for Radius in target_R:
	new_theta,new_phi,new_data=R_slice(Radius,coord_logR,coord_theta,coord_phi,data)
	new_data=data_slice(new_data)
	interp_grid=griddata((new_theta,new_phi), new_data, (grid_theta,grid_phi), method='linear')
	data_list.append(interp_grid)
	if min_data==None:
		min_data=min(interp_grid.flatten())
		max_data=max(interp_grid.flatten())
	else:
		min_data=min(min(interp_grid.flatten()),min_data)
		max_data=max(max(interp_grid.flatten()),max_data)	

print(min_data,max_data)


X_grid_1=np.sin(grid_theta)*np.cos(grid_phi)
Y_grid_1=np.sin(grid_theta)*np.sin(grid_phi)
Z_grid_1=np.cos(grid_theta)

call_result=call(["mkdir","./anim_temp"])


mlab.figure(bgcolor=(0.95,0.95,0.95),size=(1200,1000))

#Do the 0th frame twice, first one has a problem
mlab.clf()

sphere_mesh=mlab.mesh(X_grid_1[::-1,:]*target_R[0],Y_grid_1[::-1,:]*target_R[0],Z_grid_1[:,:]*target_R[0],scalars=data_list[0][:,:],colormap='hot',vmin=0.0,vmax=1.0)
sphere_mesh.actor.property.backface_culling = True
sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True
#mlab.colorbar(orientation="vertical")

mlab.view(azimuth=0, elevation=130, roll=90, distance=6.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)

mlab.savefig("./anim_temp/img{:03d}.png".format(0))

#0th frame properly
plot_idx=0
for idx_f in range(len(target_R)):
	mlab.clf()

	sphere_mesh=mlab.mesh(X_grid_1[::-1,:]*target_R[idx_f],Y_grid_1[::-1,:]*target_R[idx_f],Z_grid_1[:,:]*target_R[idx_f],scalars=data_list[idx_f][:,:],colormap='RdPu',vmin=0.0,vmax=1.0)
	sphere_mesh.actor.property.backface_culling = True
	sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = False
	mlab.colorbar(orientation="vertical")

	if plot_field_lines:
		for idx2 in range(min(len(R_start),len(theta_start),len(phi_start))):
			current_line_X=field_line_X[idx2]
			current_line_Y=field_line_Y[idx2]
			current_line_Z=field_line_Z[idx2]
			for idx3 in range(len(current_line_X)-1,-1,-1):
				if field_line_R[idx2][idx3]<target_R[idx_f]:
					current_line_X=np.delete(current_line_X,idx3)
					current_line_Y=np.delete(current_line_Y,idx3)
					current_line_Z=np.delete(current_line_Z,idx3)
			if len(current_line_X)>0:
				mlab.plot3d(current_line_X,current_line_Y,current_line_Z,line_width=0.01,color=(0,1,0),tube_radius=0.004)

	mlab.view(azimuth=0, elevation=130, roll=90, distance=6.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)

	mlab.savefig("./anim_temp/img{:03d}.png".format(plot_idx))
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(1,frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(target_R)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Flicks_R3D.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



