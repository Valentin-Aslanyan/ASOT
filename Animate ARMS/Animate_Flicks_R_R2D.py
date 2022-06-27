
"""
Animate flicks data at different radii, given time, in 2D Carrington coordinates

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""

target_R=[1.0,1.05,1.1,1.15,1.2,1.25]	
file_directory="./"
flicks_file="flicks.0000000"

frames_per_step=10
frames_per_sec=2
pad_start_frames=2
pad_end_frames=2

phi_limits=[-50,70]	#None for default
theta_limits=[10,60]	#degrees

#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return np.sqrt((data_in[:,9]**2+data_in[:,10]**2+data_in[:,11]**2)/(data_in[:,6]**2+data_in[:,7]**2+data_in[:,8]**2))

def R_label_function(R):
	return "$R={:.2f}$".format(R)


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call
from scipy.interpolate import griddata


size_t_grid=400
size_p_grid=500


if theta_limits!=None:
	grid_theta=np.linspace((theta_limits[0]+90.0)*DEG2RAD,(theta_limits[1]+90.0)*DEG2RAD,num=size_t_grid)
else:
	grid_theta=np.linspace(0.0,np.pi,num=size_t_grid)
if phi_limits!=None:
	grid_phi=np.linspace(phi_limits[0]*DEG2RAD,phi_limits[1]*DEG2RAD,num=size_p_grid)
else:
	grid_phi=np.linspace(-np.pi,np.pi,num=size_p_grid)
grid_theta,grid_phi=np.meshgrid(grid_theta,grid_phi)
data_list=[]
min_data=None
max_data=None
time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)
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


call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(target_R)):
	plt.clf()
	fig=plt.figure(figsize=(14,10))
	ax=fig.gca()
	plt.tight_layout()
	color_plot=plt.pcolormesh(grid_phi*RAD2DEG,grid_theta*RAD2DEG+90.0,data_list[idx_f],cmap='hot_r',vmin=min_data,vmax=1.0)#max_data)
	cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
	cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
	cbar.set_label(label=r"Quantity",fontsize=20)
	plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
	plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	if phi_limits!=None:
		if theta_limits!=None:
			time_label=plt.text(phi_limits[0]+0.1*(phi_limits[1]-phi_limits[0]),theta_limits[0]+0.9*(theta_limits[1]-theta_limits[0]),R_label_function(target_R[idx_f]),fontsize=26)
		else:
			time_label=plt.text(phi_limits[0]+0.1*(phi_limits[1]-phi_limits[0]),70,R_label_function(target_R[idx_f]),fontsize=26)
	elif theta_limits!=None:
		time_label=plt.text(-170,theta_limits[0]+0.9*(theta_limits[1]-theta_limits[0]),R_label_function(target_R[idx_f]),fontsize=26)			
	else:
		time_label=plt.text(-170,70,R_label_function(target_R[idx_f]),fontsize=26)
	if phi_limits!=None:
		plt.xlim(phi_limits)
	if theta_limits!=None:
		plt.ylim(theta_limits)

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
	if idx_f==len(target_R)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Flicks_R2D.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



