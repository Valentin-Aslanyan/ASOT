
"""
Display squashing factor from QSL Squasher (not included in this distribution) on 3D sphere
Optionally display fieldlines, magnetic nulls and and outer translucent sphere showing outside squashing factor
"""

file_directory="./PFLS/0000000/"
R_start    =[ +1.000, +1.000, +1.000, +1.000]
theta_start=[+24.133,+23.503,+22.930,+22.472]
phi_start  =[ -8.881, -8.938, -8.594, -7.620]
plot_R3=True
plot_nulls=True


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
if plot_nulls:
	import HQVseg
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

if plot_nulls:
	null_locations=HQVseg.nulls_from_dir(file_directory)
	null_point_X=null_locations[:,0]*np.sin(null_locations[:,1]/180.0*np.pi)*np.cos(null_locations[:,2]/180.0*np.pi)
	null_point_Y=-null_locations[:,0]*np.sin(null_locations[:,1]/180.0*np.pi)*np.sin(null_locations[:,2]/180.0*np.pi)
	null_point_Z=null_locations[:,0]*np.cos(null_locations[:,1]/180.0*np.pi)



fig=plt.figure("",figsize=(10,7))
color_plot=plt.pcolormesh(phi_grid1,theta_grid1,Q_grid1[:-1,:-1],cmap='RdBu_r',vmin=-10,vmax=10)


mlab.figure(bgcolor=(1,1,1))

sphere_mesh=mlab.mesh(X_grid1[::-1,:],Y_grid1[::-1,:],Z_grid1[:,:],scalars=Q_grid1[:,:],colormap='RdBu',vmin=-10,vmax=10)
sphere_mesh.actor.property.backface_culling = True
sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True

if plot_R3:
	sphere_mesh3=mlab.mesh(X_grid3[::-1,:],Y_grid3[::-1,:],Z_grid3[:,:],scalars=Q_grid3[:,:],colormap='RdBu',vmin=-10,vmax=10,opacity=0.4)
	sphere_mesh3.actor.property.backface_culling = True
	sphere_mesh3.module_manager.scalar_lut_manager.reverse_lut = True


for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
	th_temp,ph_temp=change_angular_coords(theta_start[idx],phi_start[idx],from_type='carrington',to_type='flicks')
	field_line_start=np.array([R_start[idx],th_temp+0.5*np.pi,ph_temp])
	field_line_sph=field_line_spherical(field_line_start,R,theta,phi,B,1.0,2.9,step_size=1E-2)
	field_line_X=field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2])
	field_line_Y=-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2])
	field_line_Z=field_line_sph[:,0]*np.cos(field_line_sph[:,1])
	mlab.plot3d(field_line_X,field_line_Y,field_line_Z,line_width=0.01,color=(0.22,0.424,0.69),tube_radius=0.004)#color=(0,0,0)
	plt.plot(field_line_sph[:,2],field_line_sph[:,1],color="black")


if plot_nulls:
	null_point=mlab.points3d(null_point_X,null_point_Y,null_point_Z,scale_factor=0.025,color=(0.0, 1.0, 0.0))
	null_point.actor.property.backface_culling = True
	plt.plot(null_locations[:,2]*DEG2RAD,null_locations[:,1]*DEG2RAD,'+',color="#9933FF")


mlab.view(azimuth=0, elevation=110, roll=90, distance=4.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)

mlab.show()
cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)
plt.show()



