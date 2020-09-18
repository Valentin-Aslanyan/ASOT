

target_R=1.1	#
file_directory="./"
file_labels=['0044295',
 '0045393',
 '0046491',
 '0047592',
 '0048697',
 '0049811',
 '0050943',
 '0052102',
 '0053303',
 '0054553',
 '0055843',
 '0057191',
 '0058693',
 '0060485',
 '0062657',
 '0064977',
 '0067537',
 '0070299',
 '0073197',
 '0076220',
 '0079415',
 '0082818',
 '0086411',
 '0090247',
 '0094203',
 '0098514',
 '0103182',
 '0108224',
 '0113626',
 '0119721',
 '0125951',
 '0132228',
 '0138526',
 '0144829',
 '0151131',
 '0157436',
 '0163745',
 '0170056',
 '0176371',
 '0182686',
 '0189003',
 '0195323',
 '0201646',
 '0207970',
 '0214295',
 '0220620',
 '0226944',
 '0233269',
 '0239596',
 '0245921',
 '0252247',
 '0258573',
 '0264899',
 '0271225',
 '0277550',
 '0283876',
 '0290202',
 '0296528',
 '0302854',
 '0309180',
 '0315507',
 '0321832',
 '0328157',
 '0334482',
 '0340807',
 '0347131',
 '0353456',
 '0359780',
 '0366104',
 '0372428',
 '0378751',
 '0385075',
 '0391400',
 '0397724',
 '0404049',
 '0410375',
 '0416699',
 '0423024',
 '0429348',
 '0435672',
 '0441995',
 '0448319',
 '0454643',
 '0460966',
 '0467288',
 '0473610',
 '0479932',
 '0486254',
 '0492576',
 '0498898',
 '0505219']

frames_per_step=2
frames_per_sec=2
pad_start_frames=2
pad_end_frames=2

phi_limits=[-50,70]	#None for default
theta_limits=[100,150]	#degrees

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


if theta_limits!=None:
	grid_theta=np.linspace(theta_limits[0]/180.0*np.pi,theta_limits[1]/180.0*np.pi,num=300)
else:
	grid_theta=np.linspace(0.0,np.pi,num=300)
if phi_limits!=None:
	grid_phi=np.linspace(phi_limits[0]/180.0*np.pi,phi_limits[1]/180.0*np.pi,num=200)
else:
	grid_phi=np.linspace(-np.pi,np.pi,num=200)
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
	plt.title("log$_{10}$ Mass density, $R=1.1$",fontsize=20)
	color_plot=plt.pcolormesh(grid_phi*180.0/np.pi,grid_theta*180.0/np.pi,data_list[idx_f],cmap='pink_r',vmin=min_data,vmax=max_data)
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
			time_label=plt.text(phi_limits[0]+0.1*(phi_limits[1]-phi_limits[0]),160,time_label_function(time_list[idx_f]),fontsize=26)
	elif theta_limits!=None:
		time_label=plt.text(-170,theta_limits[0]+0.9*(theta_limits[1]-theta_limits[0]),time_label_function(time_list[idx_f]),fontsize=26)			
	else:
		time_label=plt.text(-170,160,time_label_function(time_list[idx_f]),fontsize=26)
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



