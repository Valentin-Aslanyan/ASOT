
"""
Animate flicks data at given longitude (phi) in 2D Carrington coordinates, optionally with fieldlines

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""

target_phi=0.0	#degrees
file_directory="./"
file_labels=['0000000']

R_limits=None		#None or Solar radii

plot_gridfieldlines=False
gridfieldlines_lims=[1.0,3.0]
R_start    =[ +1.000]
theta_start=[+11.471]
phi_start  =[ +0.000]

frames_per_step=3
frames_per_sec=10
pad_start_frames=3
pad_end_frames=3

#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return np.log((data_in[:,9]**2+data_in[:,10]**2+data_in[:,11]**2)/(data_in[:,6]**2+data_in[:,7]**2+data_in[:,8]**2))


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call
from scipy.interpolate import griddata


def is_block_in_limits(coord_logR,coord_theta,coord_phi,target_phi,R_limits):
	in_limits=False
	if R_limits!=None:
		if coord_phi[0]<=target_phi*DEG2RAD and coord_phi[1]>=target_phi*DEG2RAD and np.exp(coord_logR[1])<=R_limits[1] and np.exp(coord_logR[0])>=R_limits[0]:
			in_limits=True
	else:
		if coord_phi[0]<=target_phi*DEG2RAD and coord_phi[1]>=target_phi*DEG2RAD:
			in_limits=True
	return in_limits

grid_x=[]
grid_y=[]
plot_data_grid=[]
arc_outer=[]
radial_grid=[]
arc1=[]
arc2=[]
fieldlines=[]
plot_data_max=0.0
num_arc=200

for idx_f in range(len(file_labels)):
	time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,'flicks.'+file_labels[idx_f])
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
	if plot_data_max==0.0:
		plot_data_max=max(target_data.flatten())
	else:
		plot_data_max=max(plot_data_max,max(target_data.flatten()))

	min_new_logR=min(new_logR)
	max_new_logR=max(new_logR)
	min_new_theta=min(new_theta)
	max_new_theta=max(new_theta)
	min_newR2=np.exp(min_new_logR)**2

	grid_x_curr=np.linspace(min(target_x),max(target_x),num=500)
	grid_y_curr=np.linspace(min(target_y),max(target_y),num=500)
	grid_x_curr,grid_y_curr=np.meshgrid(grid_x_curr,grid_y_curr)
	plot_data_grid_curr=griddata((target_x,target_y), target_data, (grid_x_curr,grid_y_curr), method='linear')

	for idx_x in range(len(plot_data_grid_curr[0,:])):
		for idx_y in range(len(plot_data_grid_curr[:,0])):
			if grid_x_curr[0,idx_x]**2+grid_y_curr[idx_y,0]**2<min_newR2:
				plot_data_grid_curr[idx_y,idx_x]=np.nan
	grid_x.append(grid_x_curr)
	grid_y.append(grid_y_curr)
	plot_data_grid.append(plot_data_grid_curr)

	arc_r=np.concatenate((np.exp(min_new_logR)*np.ones((num_arc)),np.exp(max_new_logR)*np.ones((num_arc)),[np.exp(min_new_logR)]))
	arc_theta=np.concatenate((np.linspace(min_new_theta,max_new_theta,num=num_arc),np.linspace(max_new_theta,min_new_theta,num=num_arc),[min_new_theta]))
	arc_outer_curr=np.zeros((2,len(arc_theta)))
	arc_outer_curr[0,:]=arc_r*np.sin(arc_theta)
	arc_outer_curr[1,:]=arc_r*np.cos(np.pi-arc_theta)
	arc_outer.append(arc_outer_curr)

	if plot_gridfieldlines:
		nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
		n1pm1=n1p-1
		n2pm1=n2p-1
		n3pm1=n3p-1
		B_flicks=data[:,:,:,:,6:9]
		B_flicks[:,:,:,:,1]*=-1.0
		num_valid_blocks=0
		for idx in range(len(coord_logR)):
			if is_block_in_limits(coord_logR[idx],coord_theta[idx],coord_phi[idx],target_phi,R_limits):
				num_valid_blocks+=1
		radial_grid_curr=np.zeros((num_valid_blocks,2,4))
		arc1_curr=np.zeros((num_valid_blocks,2,100))
		arc2_curr=np.zeros((num_valid_blocks,2,100))
		idx_b=0
		for idx in range(len(coord_logR)):
			if is_block_in_limits(coord_logR[idx],coord_theta[idx],coord_phi[idx],target_phi,R_limits):
				radial_grid_curr[idx_b,0,0]=np.exp(coord_logR[idx][0])*np.sin(coord_theta[idx][0])
				radial_grid_curr[idx_b,1,0]=np.exp(coord_logR[idx][0])*np.cos(np.pi-coord_theta[idx][0])
				radial_grid_curr[idx_b,0,1]=np.exp(coord_logR[idx][1])*np.sin(coord_theta[idx][0])
				radial_grid_curr[idx_b,1,1]=np.exp(coord_logR[idx][1])*np.cos(np.pi-coord_theta[idx][0])
				radial_grid_curr[idx_b,0,2]=np.exp(coord_logR[idx][1])*np.sin(coord_theta[idx][1])
				radial_grid_curr[idx_b,1,2]=np.exp(coord_logR[idx][1])*np.cos(np.pi-coord_theta[idx][1])
				radial_grid_curr[idx_b,0,3]=np.exp(coord_logR[idx][0])*np.sin(coord_theta[idx][1])
				radial_grid_curr[idx_b,1,3]=np.exp(coord_logR[idx][0])*np.cos(np.pi-coord_theta[idx][1])
				anglearc=np.linspace(coord_theta[idx][0],coord_theta[idx][1],num=100)
				arc1_curr[idx_b,0,:]=np.exp(coord_logR[idx][0])*np.sin(anglearc)
				arc1_curr[idx_b,1,:]=np.exp(coord_logR[idx][0])*np.cos(np.pi-anglearc)
				arc2_curr[idx_b,0,:]=np.exp(coord_logR[idx][1])*np.sin(anglearc)
				arc2_curr[idx_b,1,:]=np.exp(coord_logR[idx][1])*np.cos(np.pi-anglearc)
				idx_b+=1
		radial_grid.append(radial_grid_curr)
		arc1.append(arc1_curr)
		arc2.append(arc2_curr)
		
		fieldlines_curr=[]
		for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
			th_temp,ph_temp=change_angular_coords(theta_start[idx],phi_start[idx],from_type='carrington',to_type='flicks')
			field_line_start=np.array([R_start[idx],th_temp+0.5*np.pi,ph_temp])
			field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
			field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
			field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
			field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
			fieldlines_curr.append([field_line_X,field_line_Z])
		fieldlines.append(fieldlines_curr)

call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(file_labels)):
	plt.clf()
	fig=plt.figure("Output",figsize=(3.1*1.5,6*1.5))
	plt.pcolormesh(grid_x[idx_f],grid_y[idx_f],plot_data_grid[idx_f],cmap='RdPu',rasterized=True,vmax=plot_data_max)#,vmin=0.0,vmax=1.0)
	plt.plot(arc_outer[idx_f][0,:],arc_outer[idx_f][1,:],color="blue",linewidth=2)
	if plot_gridfieldlines:
		for idx in range(len(arc1[idx_f][:,0,0])):
			plt.plot(radial_grid[idx_f][idx,0,:2],radial_grid[idx_f][idx,1,:2],color="grey",linewidth=1)
			plt.plot(radial_grid[idx_f][idx,0,2:],radial_grid[idx_f][idx,1,2:],color="grey",linewidth=1)
			plt.plot(arc1[idx_f][idx,0,:],arc1[idx_f][idx,1,:],color="grey",linewidth=1)
			plt.plot(arc2[idx_f][idx,0,:],arc2[idx_f][idx,1,:],color="grey",linewidth=1)
		for idx in range(len(fieldlines[idx_f])):
			plt.plot(fieldlines[idx_f][idx][0],fieldlines[idx_f][idx][1],color="black",linewidth=2)

	plt.axis('off')
	plt.xlim([0,3.1])
	plt.ylim([-3,3])

	plt.savefig("./anim_temp/img{:03d}.png".format(plot_idx), format="png", dpi=100,bbox_inches='tight',pad_inches=0.1)
	plt.close(fig)
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(1,frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(file_labels)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Flicks_t_Phi2D.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



