

file_directory="./PFLS/0044295/"
R_start    =[+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000]
theta_start=[+2.510,+2.028,+2.048,+2.054,+1.982,+1.990,+1.951,+1.961,+2.006,+2.012,+1.891,+1.834,+1.798,+1.740,+1.692,+1.834,+2.021,+2.706,+2.694,+2.657,+2.543,+2.550,+2.513,+2.527]
phi_start  =[-0.811,+1.135,-0.282,-0.282,-0.015,-0.015,+0.215,+0.216,-0.162,-0.155,+0.546,+0.318,+0.160,-0.031,-0.202,-0.277,-0.311,-0.210,-0.211,+0.059,+0.320,+0.333,+0.727,+0.740]

plot_R3=True


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from mayavi import mlab


R1,theta_grid1,phi_grid1,Q1=parse_QSL_Rbinfile(file_directory+'qslR1.bin')
Q_grid1=np.sign(Q1)*np.log(abs(Q1))
Q_grid1[np.isinf(Q_grid1)]=np.nan
X_grid1=np.sin(theta_grid1)*np.cos(phi_grid1)
Y_grid1=np.sin(theta_grid1)*np.sin(phi_grid1)
Z_grid1=np.cos(theta_grid1)

if plot_R3:
	R3,theta_grid3,phi_grid3,Q3=parse_QSL_Rbinfile(file_directory+'qslR3.bin')
	Q_grid3=np.sign(Q3)*np.log(abs(Q3))
	Q_grid3[np.isinf(Q_grid3)]=np.nan
	X_grid3=2.89*np.sin(theta_grid3)*np.cos(phi_grid3)
	Y_grid3=2.89*np.sin(theta_grid3)*np.sin(phi_grid3)
	Z_grid3=2.89*np.cos(theta_grid3)


phi,theta,R,B=parse_QSL_folder(file_directory)
phi=(phi-180.0)*np.pi/180.0
theta=(theta+90.0)*np.pi/180.0


#field_line_start=np.array([1.0,1.99,-0.1])
#field_line_sph=field_line_spherical(field_line_start,R,theta,phi,B,1.0,3.0,step_size=1E-2)
#field_line_X=field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2])
#field_line_Y=-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2])
#field_line_Z=field_line_sph[:,0]*np.cos(field_line_sph[:,1])

#plt.figure()
#print(np.shape(Q_grid1),np.shape(phi_grid1),np.shape(theta_grid1))
#plt.plot(phi_grid1[:,0])
#plt.plot(theta_grid1[0,:])


fig=plt.figure("",figsize=(10,7))
color_plot=plt.pcolormesh(phi_grid1,theta_grid1,Q_grid1,cmap='RdBu_r',vmin=-10,vmax=10)


#mlab.figure(bgcolor=(1,0.8,0.55))
mlab.figure(bgcolor=(1,1,1))

sphere_mesh=mlab.mesh(X_grid1[::-1,:],Y_grid1[::-1,:],Z_grid1[:,:],scalars=Q_grid1[:,:],colormap='RdBu',vmin=-10,vmax=10)
sphere_mesh.actor.property.backface_culling = True
sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True
#mlab.colorbar(orientation="vertical")

if plot_R3:
	sphere_mesh3=mlab.mesh(X_grid3[::-1,:],Y_grid3[::-1,:],Z_grid3[:,:],scalars=Q_grid3[:,:],colormap='RdBu',vmin=-10,vmax=10,opacity=0.4)
	sphere_mesh3.actor.property.backface_culling = True
	sphere_mesh3.module_manager.scalar_lut_manager.reverse_lut = True


for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
	field_line_start=np.array([R_start[idx],theta_start[idx],phi_start[idx]])
	field_line_sph=field_line_spherical(field_line_start,R,theta,phi,B,1.0,2.9,step_size=1E-2)
	field_line_X=field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2])
	field_line_Y=-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2])
	field_line_Z=field_line_sph[:,0]*np.cos(field_line_sph[:,1])
	mlab.plot3d(field_line_X,field_line_Y,field_line_Z,line_width=0.01,color=(0.22,0.424,0.69),tube_radius=0.004)#color=(0,0,0)
	plt.plot(field_line_sph[:,2],field_line_sph[:,1],color="black")


mlab.view(azimuth=0, elevation=110, roll=90, distance=4.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)

mlab.show()
cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)
plt.show()



