
"""
Animate successive slices through flicks files with data interpolated along a given field line

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""

target_phi=0.0	#degrees
file_directory="./"
flicks_files=["flicks.0000000"]

frames_per_step=3
frames_per_sec=2
pad_start_frames=6
pad_end_frames=2

R_start    =[ +1.000]
theta_start=[+11.471]
phi_start  =[ +0.000]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call



field_lines=[]
grid_coords=[]
output_times=[]
max_J=0.0
max_JoverB=0.0
max_size_field_line=0

for fname in flicks_files:

	time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,fname)
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	B_flicks=data[:,:,:,:,6:9]
	B_flicks[:,:,:,:,1]*=-1.0
	data_flicks=np.zeros((np.shape(data)[0],np.shape(data)[1],np.shape(data)[2],np.shape(data)[3],3))
	data_flicks[:,:,:,:,0]=np.sqrt((data[:,:,:,:,9]**2+data[:,:,:,:,10]**2+data[:,:,:,:,11]**2))
	data_flicks[:,:,:,:,1]=np.sqrt((data[:,:,:,:,9]**2+data[:,:,:,:,10]**2+data[:,:,:,:,11]**2)/(data[:,:,:,:,6]**2+data[:,:,:,:,7]**2+data[:,:,:,:,8]**2))
	data_flicks[:,:,:,:,2]=np.sqrt((data[:,:,:,:,6]**2+data[:,:,:,:,7]**2+data[:,:,:,:,8]**2))

	output_times.append(time)
	grid_coords_current=[]
	for idx in range(len(coord_logR)):
		if coord_phi[idx][0]<=target_phi/180.0*np.pi and coord_phi[idx][1]>=target_phi/180.0*np.pi:
			grid_coords_current.append([coord_logR[idx],coord_theta[idx]])
	grid_coords.append(grid_coords_current)
	field_lines_current=[]
	for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
		th_temp,ph_temp=change_angular_coords(theta_start[idx],phi_start[idx],from_type='carrington',to_type='flicks')
		field_line_start=np.array([R_start[idx],th_temp+0.5*np.pi,ph_temp])
		field_line_sph=data_along_field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,data_flicks,1.0,2.9,nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		max_J=max(max_J,max(field_line_sph[:,3]))
		max_JoverB=max(max_JoverB,max(field_line_sph[:,4]))
		field_lines_current.append(field_line_sph)
		max_size_field_line=max(max_size_field_line,len(field_line_sph[:,0]))
	field_lines.append(field_lines_current)


call_result=call(["mkdir","./anim_temp"])

plot_idx=0
for idx_f in range(len(flicks_files)):
	plt.clf()
	fig1=plt.figure(figsize=(12,10))
	plt.subplot(1,2,1)
	ax1_1=fig1.gca()
	plt.title("$t={:d}$ s".format(int(output_times[idx_f])),fontsize=20)
	for idx in range(len(grid_coords[idx_f])):
		x1=np.exp(grid_coords[idx_f][idx][0][0])*np.sin(grid_coords[idx_f][idx][1][0])
		z1=np.exp(grid_coords[idx_f][idx][0][0])*np.cos(np.pi-grid_coords[idx_f][idx][1][0])
		x2=np.exp(grid_coords[idx_f][idx][0][1])*np.sin(grid_coords[idx_f][idx][1][0])
		z2=np.exp(grid_coords[idx_f][idx][0][1])*np.cos(np.pi-grid_coords[idx_f][idx][1][0])
		x3=np.exp(grid_coords[idx_f][idx][0][1])*np.sin(grid_coords[idx_f][idx][1][1])
		z3=np.exp(grid_coords[idx_f][idx][0][1])*np.cos(np.pi-grid_coords[idx_f][idx][1][1])
		x4=np.exp(grid_coords[idx_f][idx][0][0])*np.sin(grid_coords[idx_f][idx][1][1])
		z4=np.exp(grid_coords[idx_f][idx][0][0])*np.cos(np.pi-grid_coords[idx_f][idx][1][1])
		plt.plot([x1,x2],[z1,z2],color="grey",linewidth=1)
		plt.plot([x3,x4],[z3,z4],color="grey",linewidth=1)
		anglearc=np.linspace(grid_coords[idx_f][idx][1][0],grid_coords[idx_f][idx][1][1],num=100)
		arc1x=np.exp(grid_coords[idx_f][idx][0][0])*np.sin(anglearc)
		arc1y=np.exp(grid_coords[idx_f][idx][0][0])*np.cos(np.pi-anglearc)
		arc2x=np.exp(grid_coords[idx_f][idx][0][1])*np.sin(anglearc)
		arc2y=np.exp(grid_coords[idx_f][idx][0][1])*np.cos(np.pi-anglearc)
		plt.plot(arc1x,arc1y,color="grey",linewidth=1)
		plt.plot(arc2x,arc2y,color="grey",linewidth=1)
	for idx in range(len(field_lines[idx_f])):
		field_line_X=field_lines[idx_f][idx][:,0]*np.sin(np.pi-field_lines[idx_f][idx][:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_lines[idx_f][idx][:,0]*np.sin(np.pi-field_lines[idx_f][idx][:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_lines[idx_f][idx][:,0]*np.cos(np.pi-field_lines[idx_f][idx][:,1])
		plt.plot(field_line_X,field_line_Z,color="black",linewidth=2)
	plt.xlim([0,3.1])
	plt.ylim([-3,3])
	plt.xlabel(r'$R$ [$R_{\odot}$]',fontsize=20)
	plt.ylabel(r'$Z$ [$R_{\odot}$]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

	ax1_2=fig1.add_axes([0.52, 0.2, 0.26, 0.26])
	plt.title(r"Normalized $|\vec{J}|$",fontsize=20)
	plt.plot(field_lines[idx_f][0][:,3]/max_J,color="black",linewidth=1)
	plt.xlim([0,max_size_field_line])
	plt.ylim([0.0,1.05])
	plt.xlabel(r'Path along field line [a.u.]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

	ax1_3=fig1.add_axes([0.52, 0.6, 0.26, 0.26])
	plt.title(r"Normalized $|\vec{J}|/|\vec{B}|$",fontsize=20)
	plt.plot(field_lines[idx_f][0][:,4]/max_JoverB,color="black",linewidth=1)
	plt.xlim([0,max_size_field_line])
	plt.ylim([0.0,1.05])
	plt.xlabel(r'Path along field line [a.u.]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)


	plt.savefig("./anim_temp/img{:03d}.png".format(plot_idx), format="png", dpi=100,bbox_inches='tight',pad_inches=0.1)
	plt.close(fig)
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(flicks_files)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_AlongFieldline.mp4"])
call_result=call(["rm","-r","./anim_temp/"])





