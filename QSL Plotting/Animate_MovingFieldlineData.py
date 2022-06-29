

target_phi=-8.6	#degrees
solar_Radius=1.0*7E10
phi_limits=[-20,0]
theta_limits=[100,120]
target_directory="./"
plot_cage=True
use_Fortran_library=True
SurfaceLine_initial="./SurfaceLineData/SLine_4000.bin"
t_start=4.0e+3
delta_t=10.0

frames_per_step=3
frames_per_sec=3
pad_start_frames=2
pad_end_frames=2

file_IDs=["0044295",
"0046491",
"0047589",
"0048687",
"0049786",
"0050887",
"0051993",
"0053109",
"0054245",
"0055411",
"0056624",
"0057883",
"0059184",
"0060557",
"0062111",
"0063969",
"0066124",
"0068429",
"0070944",
"0073595",
"0076407",
"0079371"]

advect_startpoints=True
R_start    =[+1.000]
theta_start=[+2.010]
phi_start  =[-0.121]

fieldline_type_colors=["#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666"]

import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call


num_startpoints=min(len(R_start),len(theta_start),len(phi_start))
Time_profiles,sflow_parameters=parse_sflow_CNT(os.path.join(target_directory,"arms.cnt"))
t=t_start
cage_stpts=np.array([phi_start[:num_startpoints],theta_start[:num_startpoints]])
cage_stpts[0,:]=cage_stpts[0,:]/np.pi
cage_stpts[1,:]=1.0-cage_stpts[1,:]/np.pi
cage_stpts_output=np.zeros((2,len(cage_stpts[0,:]),len(file_IDs)))


field_lines=[]
grid_coords=[]
output_times=[]
max_J=0.0
max_JoverB=0.0
max_B=0.0
max_rho=0.0
max_v=0.0
max_JdotB=0.0
max_size_field_line=0
Q_list=[]
for idx_f in range(len(file_IDs)):
	R,theta,phi,Q=parse_QSL_Rbinfile(os.path.join(target_directory,'PFLS',file_IDs[idx_f],'qslR1.bin'))
	Q_grid=np.sign(Q)*np.log(abs(Q))
	Q_grid[np.isinf(Q_grid)]=0.0#np.nan		#Note: animate bug if nan in data
	Q_list.append(Q_grid)

	time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(target_directory,"flicks."+file_IDs[idx_f])
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	B_flicks=data[:,:,:,:,6:9]
	B_flicks[:,:,:,:,1]*=-1.0
	data_flicks=np.zeros((np.shape(data)[0],np.shape(data)[1],np.shape(data)[2],np.shape(data)[3],6))
	data_flicks[:,:,:,:,0]=np.sqrt((data[:,:,:,:,9]**2+data[:,:,:,:,10]**2+data[:,:,:,:,11]**2))
	data_flicks[:,:,:,:,1]=np.sqrt((data[:,:,:,:,9]**2+data[:,:,:,:,10]**2+data[:,:,:,:,11]**2)/(data[:,:,:,:,6]**2+data[:,:,:,:,7]**2+data[:,:,:,:,8]**2))
	data_flicks[:,:,:,:,2]=np.sqrt((data[:,:,:,:,6]**2+data[:,:,:,:,7]**2+data[:,:,:,:,8]**2))
	data_flicks[:,:,:,:,3]=data[:,:,:,:,0]
	data_flicks[:,:,:,:,4]=np.sqrt((data[:,:,:,:,1]**2+data[:,:,:,:,2]**2+data[:,:,:,:,3]**2))
	data_flicks[:,:,:,:,5]=(data[:,:,:,:,9]*data[:,:,:,:,6]+data[:,:,:,:,10]*data[:,:,:,:,7]+data[:,:,:,:,11]*data[:,:,:,:,8])/(data[:,:,:,:,6]**2+data[:,:,:,:,7]**2+data[:,:,:,:,8]**2)

	if advect_startpoints:
		num_t=int((time-t)/delta_t)
		for idx_t in range(0,num_t):
			for idx_c in range(num_startpoints):
				new_pos=surface_move_rk_4(solar_Radius,cage_stpts[1,idx_c],cage_stpts[0,idx_c],t,delta_t,Time_profiles,sflow_parameters)
				cage_stpts[0,idx_c]=new_pos[1]
				cage_stpts[1,idx_c]=new_pos[0]
			t+=delta_t
		delta_t_new=time-t
		for idx_c in range(num_startpoints):
			new_pos=surface_move_rk_4(solar_Radius,cage_stpts[1,idx_c],cage_stpts[0,idx_c],t,delta_t_new,Time_profiles,sflow_parameters)
			cage_stpts[0,idx_c]=new_pos[1]
			cage_stpts[1,idx_c]=new_pos[0]
	cage_stpts_output[:,:,idx_f]=cage_stpts[:,:]


	output_times.append(time)
	grid_coords_current=[]
	for idx in range(len(coord_logR)):
		if coord_phi[idx][0]<=target_phi/180.0*np.pi and coord_phi[idx][1]>=target_phi/180.0*np.pi:
			grid_coords_current.append([coord_logR[idx],coord_theta[idx]])
	grid_coords.append(grid_coords_current)
	field_lines_current=[]
	for idx in range(num_startpoints):
		field_line_start=np.array([R_start[idx],np.pi*(1.0-cage_stpts_output[1,idx,idx_f]),np.pi*cage_stpts_output[0,idx,idx_f]])
		field_line_sph=data_along_field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,data_flicks,1.0,3.0,nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		max_J=max(max_J,max(field_line_sph[:,3]))
		max_JoverB=max(max_JoverB,max(field_line_sph[:,4]))
		max_B=max(max_B,max(field_line_sph[:,5]))
		max_rho=max(max_rho,max(field_line_sph[:,6]))
		max_v=max(max_v,max(field_line_sph[:,7]))
		max_JdotB=max(max_JdotB,max(abs(field_line_sph[:,8])))
		field_lines_current.append(field_line_sph)
		max_size_field_line=max(max_size_field_line,len(field_line_sph[:,0]))
	field_lines.append(field_lines_current)


if plot_cage:
	if use_Fortran_library:
		SurfaceLine_base=SurfaceLine_initial.replace('.bin','')
		len_fname=len(SurfaceLine_base)
		while len_fname>0 and SurfaceLine_base[len_fname-1].isnumeric():
			len_fname-=1
		SurfaceLine_base=SurfaceLine_base[:len_fname]
		SurfaceLine_files=[]
		Fcnt_file=open("SurfaceLine_Advance.cnt","w")
		print("!arms.cnt directory",file=Fcnt_file)
		print(os.path.join(target_directory,"arms.cnt"),file=Fcnt_file)
		print("!Initial cage file",file=Fcnt_file)
		print(SurfaceLine_initial,file=Fcnt_file)
		print("!t_start",file=Fcnt_file)
		print(t_start,file=Fcnt_file)
		print("!delta_t",file=Fcnt_file)
		print(delta_t,file=Fcnt_file)
		print("!solar_Radius",file=Fcnt_file)
		print(solar_Radius,file=Fcnt_file)
		print("!number of outputs (for plot_t)",file=Fcnt_file)
		print(len(output_times),file=Fcnt_file)
		print("!plot_t, each one on new line",file=Fcnt_file)
		for time in output_times:
			print(time,file=Fcnt_file)
		print("!SurfaceLine_outputs",file=Fcnt_file)
		for time in output_times:
			SurfaceLine_files.append(SurfaceLine_base+str(int(time))+'.bin')
			print(SurfaceLine_base+str(int(time))+'.bin',file=Fcnt_file)
		Fcnt_file.close()
		call_result=call(["./SurfaceLine_Advance"])
		static_cages,dynamic_cages_temp=surfaceline_read_binary_cage(SurfaceLine_files[0])
		dynamic_cages=[]
		for idx in range(len(dynamic_cages_temp)):
			dynamic_cages.append(np.zeros((np.shape(dynamic_cages_temp[idx])[0],np.shape(dynamic_cages_temp[idx])[1],len(SurfaceLine_files))))
			dynamic_cages[idx][:,:,0]=dynamic_cages_temp[idx][:,:]
		for idx in range(1,len(SurfaceLine_files)):
			static_cages,dynamic_cages_temp=surfaceline_read_binary_cage(SurfaceLine_files[idx])
			for idx2 in range(len(dynamic_cages_temp)):
				dynamic_cages[idx2][:,:,idx]=dynamic_cages_temp[idx2][:,:]
	else:
		static_cages,dynamic_cages_temp=surfaceline_read_binary_cage(SurfaceLine_initial)
		dynamic_cages=[]
		for idx in range(len(dynamic_cages_temp)):
			dynamic_cages.append(surfaceline_advance_cage(dynamic_cages_temp[idx],t_start,output_times,delta_t,os.path.join(target_directory,"arms.cnt"),solar_Radius))



call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(file_IDs)):
	fig1=plt.figure("A",figsize=(22,18))
	plt.clf()

	ax1_10=fig1.add_axes([0.02, 0.02, 0.30, 0.3])
	plt.xlim(phi_limits)
	plt.ylim(theta_limits)
	plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
	plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)

	color_plot=plt.pcolormesh(phi*180.0/np.pi,theta*180.0/np.pi,Q_list[idx_f],cmap='RdBu_r',vmin=-10,vmax=10,rasterized=True)
	if plot_cage:
		for idx_c in range(len(dynamic_cages)):
			plt.plot(dynamic_cages[idx_c][0,:,idx_f],dynamic_cages[idx_c][1,:,idx_f],'-',color="black")
		for idx_c in range(len(static_cages)):
			plt.plot(static_cages[idx_c][0,:],static_cages[idx_c][1,:],color="black")
	plt.plot(cage_stpts_output[0,:,idx_f]*180.0,(1.0-cage_stpts_output[1,:,idx_f])*180.0,'+',ms=14,color="black")

	ax1_1=fig1.add_axes([0.02, 0.45, 0.25, 0.55])
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

	ax1_2=fig1.add_axes([0.36, 0.8, 0.2, 0.16])
	plt.title(r"Normalized $|\vec{B}|$",fontsize=20)
	plt.plot(field_lines[idx_f][0][:,5]/max_B,color="black",linewidth=1)
	plt.xlim([0,max_size_field_line])
	plt.ylim([0.0,1.05])
	plt.xlabel(r'Path along field line [a.u.]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

	ax1_3=fig1.add_axes([0.36, 0.55, 0.2, 0.16])
	plt.title(r"Normalized $|\vec{J}|$",fontsize=20)
	plt.plot(field_lines[idx_f][0][:,3]/max_J,color="black",linewidth=1)
	plt.xlim([0,max_size_field_line])
	plt.ylim([0.0,1.05])
	plt.xlabel(r'Path along field line [a.u.]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

	ax1_4=fig1.add_axes([0.36, 0.3, 0.2, 0.16])
	plt.title(r"Normalized $|\vec{J}|/|\vec{B}|$",fontsize=20)
	plt.plot(field_lines[idx_f][0][:,4]/max_JoverB,color="black",linewidth=1)
	plt.xlim([0,max_size_field_line])
	plt.ylim([0.0,1.05])
	plt.xlabel(r'Path along field line [a.u.]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

	ax1_5=fig1.add_axes([0.65, 0.8, 0.2, 0.16])
	plt.title(r"Normalized $\rho$",fontsize=20)
	plt.plot(field_lines[idx_f][0][:,6]/max_rho,color="black",linewidth=1)
	plt.xlim([0,max_size_field_line])
	plt.ylim([0.0,1.05])
	plt.xlabel(r'Path along field line [a.u.]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

	ax1_6=fig1.add_axes([0.65, 0.55, 0.2, 0.16])
	plt.title(r"Normalized $|\vec{v}|$",fontsize=20)
	plt.plot(field_lines[idx_f][0][:,7]/max_v,color="black",linewidth=1)
	plt.xlim([0,max_size_field_line])
	plt.ylim([0.0,1.05])
	plt.xlabel(r'Path along field line [a.u.]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

	ax1_7=fig1.add_axes([0.65, 0.3, 0.2, 0.16])
	plt.title(r"Normalized $\vec{J}.\vec{B}/|\vec{B}|^2$",fontsize=20)
	plt.plot(field_lines[idx_f][0][:,8]/max_JdotB,color="black",linewidth=1)
	plt.xlim([0,max_size_field_line])
	plt.ylim([-1.05,1.05])
	plt.xlabel(r'Path along field line [a.u.]',fontsize=20)
	plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)


	plt.savefig("./anim_temp/img{:03d}.png".format(plot_idx), format="png", dpi=100,bbox_inches='tight',pad_inches=0.1)
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(1,frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(file_IDs)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_AlongFieldline.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



