
"""
Animate flicks data at given radius in 2D Carrington coordinates

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""

target_R=1.1	#Solar radii
file_directory="./"
file_labels=['0000000',
 '0000001',
 '0000002',
 '0000003',
 '0000004',
 '0000005',
 '0000006']

frames_per_step=2
frames_per_sec=2
pad_start_frames=2
pad_end_frames=2

phi_limits=[-50,70]	#None for default
theta_limits=[10,60]	#degrees

#Use this to define 2D data array desired: log(rho), B^2 etc...
def data_slice(data_in):
	return np.log(data_in[:,0])/np.log(10.0)

def time_label_function(time_float):
	return "$t={:.2f}T$".format((time_float-4000)/4000.0)


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
time_list=[]
min_data=None
max_data=None
for fname in file_labels:
	time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,'flicks.'+fname)
	print(fname,time)
	new_theta,new_phi,new_data=R_slice(target_R,coord_logR,coord_theta,coord_phi,data)
	new_data=data_slice(new_data)
	interp_grid=griddata((new_theta,new_phi), new_data, (grid_theta,grid_phi), method='linear')
	data_list.append(interp_grid)
	time_list.append(time)
	if min_data==None:
		min_data=min(interp_grid.flatten())
		max_data=max(interp_grid.flatten())
	else:
		min_data=min(min(interp_grid.flatten()),min_data)
		max_data=max(max(interp_grid.flatten()),max_data)

call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(file_labels)):
	plt.clf()
	fig=plt.figure(figsize=(14,10))
	ax=fig.gca()
	plt.tight_layout()
	plt.title("$R="+"{:.1f}".format(target_R)+"$",fontsize=20)
	color_plot=plt.pcolormesh(grid_phi*RAD2DEG,grid_theta*RAD2DEG+90.0,data_list[idx_f],cmap='pink_r',vmin=min_data,vmax=max_data)
	cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
	cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
	cbar.set_label(label=r"$\log(\rho)$",fontsize=20)
	plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
	plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	if phi_limits!=None:
		if theta_limits!=None:
			time_label=plt.text(phi_limits[0]+0.1*(phi_limits[1]-phi_limits[0]),theta_limits[0]+0.9*(theta_limits[1]-theta_limits[0]),time_label_function(time_list[idx_f]),fontsize=26)
		else:
			time_label=plt.text(phi_limits[0]+0.1*(phi_limits[1]-phi_limits[0]),70,time_label_function(time_list[idx_f]),fontsize=26)
	elif theta_limits!=None:
		time_label=plt.text(-170,theta_limits[0]+0.9*(theta_limits[1]-theta_limits[0]),time_label_function(time_list[idx_f]),fontsize=26)			
	else:
		time_label=plt.text(-170,70,time_label_function(time_list[idx_f]),fontsize=26)
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
	if idx_f==len(filenames)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Flicks.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



