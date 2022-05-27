
"""
Animate successive slices through bfield files at successive phi as specified by target_phi

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""


target_phi=[-20,-18,-16,-14,-12,-10,-8,-6,-4,-2,0]	#degrees
bfield_directory="./"
bfield_file="bfield.0000000"

R_limits=[1.0,6.0]		#None or Solar radii

frames_per_step=3
frames_per_sec=10
pad_start_frames=3
pad_end_frames=3


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call
from scipy.interpolate import griddata


num_arc=200
size_x_grid=2000
size_y_grid=2000

def is_block_in_limits(coord_logR,coord_theta,coord_phi,target_phi,R_limits):
	in_limits=False
	if R_limits!=None:
		if coord_phi[0]<=target_phi/180.0*np.pi and coord_phi[1]>=target_phi/180.0*np.pi and np.exp(coord_logR[1])<=R_limits[1] and np.exp(coord_logR[0])>=R_limits[0]:
			in_limits=True
	else:
		if coord_phi[0]<=target_phi/180.0*np.pi and coord_phi[1]>=target_phi/180.0*np.pi:
			in_limits=True
	return in_limits

grid_x=[]
grid_y=[]
plot_data_grid=[]
plot_data_max=0.0
plot_data_abs_max=0.0

time,ntblks,nlblks,coord_logR,coord_theta,coord_phi,B_read=read_bfield_file(bfield_directory,bfield_file)

if R_limits==None:
	R_lower=np.exp(min(coord_logR.flatten()))
	R_upper=np.exp(max(coord_logR.flatten()))
else:
	R_lower=R_limits[0]
	R_upper=R_limits[1]
theta_lower=min(coord_theta.flatten())
theta_upper=max(coord_theta.flatten())

arc_r=np.concatenate((R_lower*np.ones((num_arc)),R_upper*np.ones((num_arc)),[R_lower]))
arc_theta=np.concatenate((np.linspace(theta_lower,theta_upper,num=num_arc),np.linspace(theta_upper,theta_lower,num=num_arc),[theta_lower]))
arc_outer=np.zeros((2,len(arc_theta)))
arc_outer[0,:]=arc_r*np.sin(arc_theta)
arc_outer[1,:]=arc_r*np.cos(np.pi-arc_theta)


for idx_f in range(len(target_phi)):
	new_logR,new_theta,B_data_in=phi_slice(target_phi[idx_f],coord_logR,coord_theta,coord_phi,B_read)
	target_x=np.exp(new_logR)*np.sin(new_theta)
	target_y=np.exp(new_logR)*np.cos(np.pi-new_theta)
	target_data=B_data_in[:,1]

	if plot_data_max==0.0:
		plot_data_max=max(target_data.flatten())
	else:
		plot_data_max=max(plot_data_max,max(target_data.flatten()))
	if plot_data_abs_max==0.0:
		plot_data_abs_max=max(abs(target_data.flatten()))
	else:
		plot_data_abs_max=max(plot_data_abs_max,max(abs(target_data.flatten())))

	grid_x_curr=np.linspace(min(target_x),max(target_x),num=size_x_grid)
	grid_y_curr=np.linspace(min(target_y),max(target_y),num=size_y_grid)
	grid_x_curr,grid_y_curr=np.meshgrid(grid_x_curr,grid_y_curr)
	plot_data_grid_curr=griddata((target_x,target_y), target_data, (grid_x_curr,grid_y_curr), method='linear')

	for idx_x in range(len(plot_data_grid_curr[0,:])):
		for idx_y in range(len(plot_data_grid_curr[:,0])):
			if grid_x_curr[0,idx_x]**2+grid_y_curr[idx_y,0]**2<R_lower**2 or grid_x_curr[0,idx_x]**2+grid_y_curr[idx_y,0]**2>R_upper**2:
				plot_data_grid_curr[idx_y,idx_x]=np.nan
	grid_x.append(grid_x_curr)
	grid_y.append(grid_y_curr)
	plot_data_grid.append(plot_data_grid_curr)


call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(target_phi)):
	plt.clf()
	fig=plt.figure("Output",figsize=(3.1*1.5,6*1.5))
	plt.title(r"$\phi="+"{:+.1f}".format(target_phi[idx_f])+"^{\circ}$",fontsize=20)
	plt.pcolormesh(grid_x[idx_f],grid_y[idx_f],plot_data_grid[idx_f],cmap='RdPu',rasterized=True,vmin=-plot_data_abs_max,vmax=plot_data_abs_max)
	plt.plot(arc_outer[0,:],arc_outer[1,:],color="blue",linewidth=2)

	plt.axis('off')
	plt.xlim([0,R_upper*1.01])
	plt.ylim([-R_upper,R_upper])

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
	if idx_f==len(target_phi)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Bfield_phi_Phi2D.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



