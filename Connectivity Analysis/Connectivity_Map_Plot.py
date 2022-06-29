
"""
Plot the squashing factor (from QSL which is not included in this repository) at the photosphere and higher up at a late time in a simulation
The connectivity map is also plotted; a selection of field lines of a given type (e.g. reconnected open) are optionally displayed at the outer boundary
"""

end_file="./PFLS/0001000/qslR1.bin"
end_file_outer="./PFLS/0001000/qslR3.bin"
connection_filename="./PFLS/Connection_0000000_0001000.bin"
start_time=0.0
end_time=1000.0
phi_limits=[-12,-2]		#None for default
theta_limits=[20,30]		#
phi_limits_out=[-60,60]		#None for default
theta_limits_out=[-10,60]	#

SurfaceLine_file="./SLine_1000.bin"

plot_top_connections=True	#Plot where the given classification(s) of field lines reach upper boundary
top_connection_types=[2]	#designations of the types of field lines to trace
top_connection_num=5	#Number of field lines to trace, or "All"
radius_inner=1.0
radius_outer=2.99


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import random


connection_colors=[(171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)]
connection_cmap=matplotlib.colors.ListedColormap(((171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)))


theta_grid,phi_grid,connection_map=load_connection_map(connection_filename)
connection_map=np.floor(connection_map/2)
R_end,theta_end,phi_end,Q_end=parse_QSL_Rbinfile(end_file)
R_out,theta_out,phi_out,Q_out=parse_QSL_Rbinfile(end_file_outer)

static_cages,dynamic_cages=surfaceline_read_binary_cage(SurfaceLine_file)


if plot_top_connections:
	phi,theta,R,B=parse_QSL_folder(os.path.dirname(end_file))
	phi=(phi-180.0)*np.pi/180.0
	theta=(theta+90.0)*np.pi/180.0
	top_connections=[]
	do_all=False
	if type(top_connection_num)==str:
		if top_connection_num.lower()=='all':
			do_all=True
	elif type(top_connection_num)!=int:
		print("Set top_connection_num to integer or 'all'")
	for idx in range(len(top_connection_types)):
		top_connections.append([])
		num_current=0
		pixel_idx=np.zeros((2,len(phi_grid[:,0])*len(theta_grid[0,:])),dtype='int')
		num_type=0
		for idx_p in range(len(phi_grid[:,0])):
			for idx_t in range(len(theta_grid[0,:])):
				if connection_map[idx_p,idx_t]==top_connection_types[idx]:
					pixel_idx[0,num_type]=idx_p
					pixel_idx[1,num_type]=idx_t
					num_type+=1
		pixel_idx=pixel_idx[:,:num_type]
		print("Connection type "+str(top_connection_types[idx])+" has "+str(num_type)+" pixels")
		if do_all:
			for idx2 in range(num_type):
				idx_p=pixel_idx[0,idx2]
				idx_t=pixel_idx[1,idx2]
				field_line_start=np.array([radius_inner,theta_grid[0,idx_t],phi_grid[idx_p,0]])
				field_line_sph=field_line_spherical(field_line_start,R,theta,phi,B,radius_inner,radius_outer,step_size=1E-2)
				if field_line_sph[-1,0]>0.9*radius_outer and field_line_sph[-1,0]<1.1*radius_outer:
					top_connections[idx].append([field_line_sph[-1,2],field_line_sph[-1,1]])
		else:
			random.seed(a=42,version=2)
			for idx2 in random.choices(list(range(num_type)),k=min(top_connection_num,num_type)):
				idx_p=pixel_idx[0,idx2]
				idx_t=pixel_idx[1,idx2]
				field_line_start=np.array([radius_inner,theta_grid[0,idx_t],phi_grid[idx_p,0]])
				field_line_sph=field_line_spherical(field_line_start,R,theta,phi,B,radius_inner,radius_outer,step_size=1E-2)
				if field_line_sph[-1,0]>0.9*radius_outer and field_line_sph[-1,0]<1.1*radius_outer:
					top_connections[idx].append([field_line_sph[-1,2],field_line_sph[-1,1]])


theta_spacing=theta_grid[0,1]-theta_grid[0,0]
phi_spacing=phi_grid[1,0]-phi_grid[0,0]
theta_spacing_out=theta_out[0,1]-theta_out[0,0]
phi_spacing_out=phi_out[1,0]-phi_out[0,0]


fig1=plt.figure("Connection map",figsize=(15,14))
plt.subplot(2,2,1)
ax1_1=fig1.gca()
plt.title(r"$R=R_{\odot}$",fontsize=20)
color_plot1_1=plt.pcolormesh((phi_grid-0.5*phi_spacing)*RAD2DEG,(theta_grid-0.5*theta_spacing)*RAD2DEG-90.0,(np.sign(Q_end)*np.log(abs(Q_end))/np.log(10.0))[:-1,:-1],cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
for idx_c in range(len(dynamic_cages)):
	plt.plot(dynamic_cages[idx_c][0,:],dynamic_cages[idx_c][1,:]-90.0,'-',color="black")
for idx_c in range(len(static_cages)):
	plt.plot(static_cages[idx_c][0,:],static_cages[idx_c][1,:]-90.0,color="black")
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)

if phi_limits!=None:
	plt.xlim(phi_limits)
if theta_limits!=None:
	plt.ylim(theta_limits)

plt.subplot(2,2,2)
ax1_2=fig1.gca()
plt.title(r"$R=R_{\mathrm{SS}}=3R_{\odot}$",fontsize=20)
color_plot1_2=plt.pcolormesh((phi_out-0.5*phi_spacing_out)*RAD2DEG,(theta_out-0.5*theta_spacing_out)*RAD2DEG-90.0,(np.sign(Q_out)*np.log(abs(Q_out))/np.log(10.0))[:-1,:-1],cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar1_2=fig1.colorbar(color_plot1_2)#,ticks=[-4,-3,-2,-1,0])
cbar1_2.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar1_2.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)
if phi_limits_out!=None:
	plt.xlim(phi_limits_out)
if theta_limits_out!=None:
	plt.ylim(theta_limits_out)

if plot_top_connections:
	for idx in range(len(top_connection_types)):
		for idx2 in range(len(top_connections[idx])):
			plt.plot(top_connections[idx][idx2][0]*RAD2DEG,top_connections[idx][idx2][1]*RAD2DEG-90.0,'H',ms=6,mec='black',mew=0.5,fillstyle='none')

plt.subplot(2,2,3)
ax1_3=fig1.gca()
color_plot1_3=plt.pcolormesh((phi_grid-0.5*phi_spacing)*RAD2DEG,(theta_grid-0.5*theta_spacing)*RAD2DEG-90.0,connection_map[:-1,:-1]+0.5,cmap=connection_cmap,vmin=0,vmax=5,rasterized=True)
for idx_c in range(len(dynamic_cages)):
	plt.plot(dynamic_cages[idx_c][0,:],dynamic_cages[idx_c][1,:]-90.0,'-',color="black")
for idx_c in range(len(static_cages)):
	plt.plot(static_cages[idx_c][0,:],static_cages[idx_c][1,:]-90.0,color="black")
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
if phi_limits!=None:
	plt.xlim(phi_limits)
if theta_limits!=None:
	plt.ylim(theta_limits)

#Key
ax1_4=fig1.add_axes([0.61, 0.2, 0.15, 0.2])
color_plot1_4=plt.pcolormesh(np.array([0,1,2]),np.array([0,1,2]),np.array([[0.5,1.5],[2.5,3.5]]),cmap=connection_cmap,vmin=0,vmax=5,rasterized=True)
ax1_4.set_xticks([])
ax1_4.set_yticks([])
plt.text(-1.3,0.55,"Connectivity",fontsize=20)
plt.text(-1.3,0.25,"retained",fontsize=20)
plt.text(-1.3,1.4,"Reconnected",fontsize=20)

plt.text(0.05,2.15,r"Open",fontsize=20)
plt.text(1.1,2.15,r"Closed",fontsize=20)
plt.text(-0.05,2.4,r"Field line end state",fontsize=20)

"""
#Old colorscheme with 8 levels
ax1_4=fig1.add_axes([0.61, 0.2, 0.26, 0.13])
color_plot1_4=plt.pcolormesh(np.array([0,1,2,3,4]),np.array([0,1,2]),np.array([[0.5,1.5,2.5,3.5,3.5],[4.5,5.5,6.5,7.5,7.5],[4.5,5.5,6.5,7.5,7.5]]),cmap=connection_cmap,vmin=0,vmax=8,rasterized=True)
ax1_4.set_xticks([])
ax1_4.set_yticks([])
plt.text(-1.3,0.4,"Convected",fontsize=20)
plt.text(-1.3,1.6,"Non-",fontsize=20)
plt.text(-1.3,1.2,"convected",fontsize=20)

plt.text(0.2,2.2,r"O$\rightarrow$O",fontsize=20)
plt.text(1.2,2.2,r"C$\rightarrow$O",fontsize=20)
plt.text(2.2,2.2,r"O$\rightarrow$C",fontsize=20)
plt.text(3.2,2.2,r"C$\rightarrow$C",fontsize=20)
plt.text(0.6,2.7,r"Open/Closed field lines",fontsize=20)
"""

plt.savefig("Connectivity_Map.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.show()



