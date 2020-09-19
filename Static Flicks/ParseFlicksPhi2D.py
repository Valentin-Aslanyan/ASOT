

target_phi=0.0	#degrees
file_directory="./"
flicks_file="flicks.0016387"
theta_lineout=None	#None or number in degrees
R_limits=None		#None or Solar radii

plot_gridfieldlines=True
gridfieldlines_lims=[1.0,3.0]
R_start4    =[+1.000,+1.000,+1.000,+1.000]
theta_start4=[+1.771]
phi_start4  =[+0.017]

R_start3    =[+1.000]
theta_start3=[+1.710]
phi_start3  =[+0.000]

R_start2    =[+1.000]
theta_start2=[+1.766]
phi_start2  =[+0.000]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)

print('Simulation time: ',time)
print(ntblks,nlblks,newgrd)
print()
print("         Min  |  Max    grid points")
print("r:     "+"{:+.3f}".format(min(coord_logR.flatten()))+" | {:+.3f}".format(max(coord_logR.flatten())))
print("theta: "+"{:+.3f}".format(min(coord_theta.flatten()))+" | {:+.3f}".format(max(coord_theta.flatten())))
print("phi:   "+"{:+.3f}".format(min(coord_phi.flatten()))+" | {:+.3f}".format(max(coord_phi.flatten())))
print()

new_logR,new_theta,new_data=phi_slice(target_phi,coord_logR,coord_theta,coord_phi,data)

def is_block_in_limits(coord_logR,coord_theta,coord_phi,target_phi,R_limits):
	in_limits=False
	if R_limits!=None:
		if coord_phi[0]<=target_phi/180.0*np.pi and coord_phi[1]>=target_phi/180.0*np.pi and np.exp(coord_logR[1])<=R_limits[1] and np.exp(coord_logR[0])>=R_limits[0]:
			in_limits=True
	else:
		if coord_phi[0]<=target_phi/180.0*np.pi and coord_phi[1]>=target_phi/180.0*np.pi:
			in_limits=True
	return in_limits
		

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
	for idx in range(min(len(R_start4),len(theta_start4),len(phi_start4))):
		field_line_start=np.array([R_start4[idx],theta_start4[idx],phi_start4[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
		plt.plot(field_line_X,field_line_Z,color="#386CB0",linewidth=2)
	for idx in range(min(len(R_start3),len(theta_start3),len(phi_start3))):
		field_line_start=np.array([R_start3[idx],theta_start3[idx],phi_start3[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
		plt.plot(field_line_X,field_line_Z,color="#FFFF99",linewidth=2)
	for idx in range(min(len(R_start2),len(theta_start2),len(phi_start2))):
		field_line_start=np.array([R_start2[idx],theta_start2[idx],phi_start2[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
		plt.plot(field_line_X,field_line_Z,color="#FDC086",linewidth=2)
	plt.axis('off')
	plt.xlim([0,3.1])
	plt.ylim([-3,3])
	plt.savefig("Grid_Field_XZ.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.figure("pflicks grid",figsize=(3.1*1.5,6*1.5))
plt.title("Grid points",fontsize=20)
plt.plot(np.exp(new_logR[:])*np.sin(new_theta[:]),np.exp(new_logR[:])*np.cos(np.pi-new_theta[:]),'o',color="red",ms=1)
plt.axis('off')
plt.xlim([0,3.1])
plt.ylim([-3,3])
plt.savefig("GridPoints_XZ.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.figure("pflicks grid logR,Theta",figsize=(7,7))
plt.title(r"Grid points $\log(r),\theta$",fontsize=20)
plt.plot(new_logR,new_theta,'o',color="red",ms=1)
plt.axis('off')

#plt.figure("pflicks grid logR,Phi",figsize=(7,7))
#plt.title(r"Grid points $\log(r),\phi$",fontsize=20)
#plt.plot(new_logR,new_phi,'o',color="red",ms=1)
#plt.axis('off')


target_x=np.exp(new_logR)*np.sin(new_theta)
target_y=np.exp(new_logR)*np.cos(np.pi-new_theta)
new_Bsquared=new_data[:,6]**2+new_data[:,7]**2+new_data[:,8]**2
new_J_B=np.sqrt((new_data[:,9]**2+new_data[:,10]**2+new_data[:,11]**2)/(new_data[:,6]**2+new_data[:,7]**2+new_data[:,8]**2))

min_new_logR=min(new_logR)
max_new_logR=max(new_logR)
min_new_theta=min(new_theta)
max_new_theta=max(new_theta)
min_newR2=np.exp(min_new_logR)**2

from scipy.interpolate import griddata
grid_x=np.linspace(min(target_x),max(target_x),num=500)
grid_y=np.linspace(min(target_y),max(target_y),num=500)
grid_x,grid_y=np.meshgrid(grid_x,grid_y)
rho_grid=griddata((target_x,target_y), np.log(new_data[:,0]), (grid_x,grid_y), method='linear')
B_grid=griddata((target_x,target_y), np.log(new_Bsquared), (grid_x,grid_y), method='linear')
JB_grid=griddata((target_x,target_y), new_J_B, (grid_x,grid_y), method='linear')
print(min(JB_grid.flatten()),max(JB_grid.flatten()))
for idx_x in range(len(B_grid[0,:])):
	for idx_y in range(len(B_grid[:,0])):
		if grid_x[0,idx_x]**2+grid_y[idx_y,0]**2<min_newR2:
			B_grid[idx_y,idx_x]=np.nan
			rho_grid[idx_y,idx_x]=np.nan
			JB_grid[idx_y,idx_x]=np.nan
num_arc=200
arc_r=np.concatenate((np.exp(min_new_logR)*np.ones((num_arc)),np.exp(max_new_logR)*np.ones((num_arc)),[np.exp(min_new_logR)]))
arc_theta=np.concatenate((np.linspace(min_new_theta,max_new_theta,num=num_arc),np.linspace(max_new_theta,min_new_theta,num=num_arc),[min_new_theta]))
arc_x=arc_r*np.sin(arc_theta)
arc_y=arc_r*np.cos(np.pi-arc_theta)


plt.figure("pcolormesh B",figsize=(3.1*1.5,6*1.5))
plt.title("$B$-field",fontsize=20)
plt.pcolormesh(grid_x,grid_y,B_grid,cmap='bone_r',rasterized=True)
plt.plot(arc_x,arc_y,color="blue",linewidth=2)
plt.axis('off')
plt.xlim([0,3.1])
plt.ylim([-3,3])
plt.savefig("B.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

if theta_lineout!=None:
	R_lineout=np.linspace(np.exp(min(new_logR)),np.exp(max(new_logR)),num=100)
	logB_lineout=griddata((new_logR,new_theta), np.log(new_Bsquared), (np.log(R_lineout),theta_lineout/180.0*np.pi*np.ones((len(R_lineout)))), method='linear')
	plt.figure("lineout B",figsize=(6.0,6.0))
	plt.plot(R_lineout,logB_lineout)


plt.figure("pcolormesh rho",figsize=(3.1*1.5,6*1.5))
plt.title("Mass density",fontsize=20)
plt.pcolormesh(grid_x,grid_y,rho_grid,cmap='hot_r',rasterized=True)
plt.plot(arc_x,arc_y,color="blue",linewidth=2)
plt.axis('off')
plt.xlim([0,3.1])
plt.ylim([-3,3])
plt.savefig("Rho.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

if theta_lineout!=None:
	R_lineout=np.linspace(np.exp(min(new_logR)),np.exp(max(new_logR)),num=100)
	logrho_lineout=griddata((new_logR,new_theta), np.log(new_data[:,0]), (np.log(R_lineout),theta_lineout/180.0*np.pi*np.ones((len(R_lineout)))), method='linear')
	plt.figure("lineout rho",figsize=(6.0,6.0))
	plt.plot(R_lineout,logrho_lineout/np.log(10.0))


plt.figure("pcolormesh J/B",figsize=(3.1*1.5,6*1.5))
plt.title(r"$\log(|\vec{J}|/|\vec{B}|)$",fontsize=20)
plt.pcolormesh(grid_x,grid_y,np.log(JB_grid),cmap='RdPu',rasterized=True)#,vmin=0.0,vmax=1.0)
plt.plot(arc_x,arc_y,color="blue",linewidth=2)
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
	for idx in range(min(len(R_start4),len(theta_start4),len(phi_start4))):
		field_line_start=np.array([R_start4[idx],theta_start4[idx],phi_start4[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
		plt.plot(field_line_X,field_line_Z,color="#386CB0",linewidth=2)
	for idx in range(min(len(R_start3),len(theta_start3),len(phi_start3))):
		field_line_start=np.array([R_start3[idx],theta_start3[idx],phi_start3[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
		plt.plot(field_line_X,field_line_Z,color="#FFFF99",linewidth=2)
	for idx in range(min(len(R_start2),len(theta_start2),len(phi_start2))):
		field_line_start=np.array([R_start2[idx],theta_start2[idx],phi_start2[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,gridfieldlines_lims[0],gridfieldlines_lims[1],nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.cos(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Y=-field_line_sph[:,0]*np.sin(np.pi-field_line_sph[:,1])#*np.sin(field_line_sph[:,2]-target_phi/180.0*np.pi)
		field_line_Z=field_line_sph[:,0]*np.cos(np.pi-field_line_sph[:,1])
		plt.plot(field_line_X,field_line_Z,color="#FDC086",linewidth=2)
plt.axis('off')
plt.xlim([0,3.1])
plt.ylim([-3,3])
plt.savefig("J_B.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


plt.show()



