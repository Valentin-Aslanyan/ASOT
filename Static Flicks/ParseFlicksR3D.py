

target_R=1.0	#degrees
file_directory="./"
flicks_file="flicks.0044295"

plot_field_lines=True  
R_start    =[+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000,+1.000]
theta_start=[+2.510,+2.028,+2.064,+1.990,+1.961,+2.012,+1.927,+1.861,+2.181,+2.659,+2.689,+2.563,+2.532,+2.514,+2.427,+2.269,+2.002]
phi_start  =[-0.811,+1.135,-0.282,-0.015,+0.216,-0.155,+0.425,+0.771,-0.517,-0.660,-0.153,+0.255,+0.401,+0.697,+0.968,+1.144,-0.162]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from mayavi import mlab


def v_Alfven(B,rho,mu=1.0):
	return 2.18E11*B/np.sqrt(rho*mu)
m_H=1.673E-24	#g



time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)
if plot_field_lines:
	nlblk,n1p,n2p,n3p,nvar=get_flicks_grid_dimensions(data)
	n1pm1=n1p-1
	n2pm1=n2p-1
	n3pm1=n3p-1
	B_flicks=data[:,:,:,:,6:9]
	B_flicks[:,:,:,:,1]*=-1.0

print('Simulation time: ',time)
print(ntblks,nlblks,newgrd)
print()
print("         Min  |  Max    grid points")
print("r:     "+"{:+.3f}".format(min(coord_logR.flatten()))+" | {:+.3f}".format(max(coord_logR.flatten())))
print("theta: "+"{:+.3f}".format(min(coord_theta.flatten()))+" | {:+.3f}".format(max(coord_theta.flatten())))
print("phi:   "+"{:+.3f}".format(min(coord_phi.flatten()))+" | {:+.3f}".format(max(coord_phi.flatten())))
print()

new_theta,new_phi,new_data=R_slice(target_R,coord_logR,coord_theta,coord_phi,data)

new_rho=new_data[:,0]
new_B=np.sqrt(new_data[:,6]**2+new_data[:,7]**2+new_data[:,8]**2)
new_v=np.sqrt(new_data[:,1]**2+new_data[:,2]**2+new_data[:,3]**2)
new_vtheta=new_data[:,2]
new_vA=new_data[:,13]

min_new_theta=min(new_theta)
max_new_theta=max(new_theta)
min_new_phi=min(new_phi)
max_new_phi=max(new_phi)

from scipy.interpolate import griddata
grid_theta=np.linspace(min_new_theta,max_new_theta,num=500)
grid_phi=np.linspace(min_new_phi,max_new_phi,num=400)
grid_theta,grid_phi=np.meshgrid(grid_theta,grid_phi)
grid_theta2=np.linspace(110/180.0*np.pi,145/180.0*np.pi,num=500)
grid_phi2=np.linspace(-40/180.0*np.pi,60/180.0*np.pi,num=400)
grid_theta2,grid_phi2=np.meshgrid(grid_theta2,grid_phi2)
v_grid=griddata((new_theta,new_phi), new_v, (grid_theta2,grid_phi2), method='linear')
B_grid=griddata((new_theta,new_phi), new_B, (grid_theta,grid_phi), method='linear')
vA_grid=griddata((new_theta,new_phi), new_vA, (grid_theta,grid_phi), method='linear')

X_grid=target_R*np.sin(grid_theta)*np.cos(grid_phi)
Y_grid=target_R*np.sin(grid_theta)*np.sin(grid_phi)
Z_grid=target_R*np.cos(grid_theta)

mlab.figure("|B| field")
sphere_mesh=mlab.mesh(X_grid[::-1,:],Y_grid[::-1,:],Z_grid,scalars=np.log(B_grid[:,:]),colormap='bone')
sphere_mesh.actor.property.backface_culling = True
sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True

if plot_field_lines:
	for idx in range(min(len(R_start),len(theta_start),len(phi_start))):
		field_line_start=np.array([R_start[idx],theta_start[idx],phi_start[idx]])
		field_line_sph=field_line_flicks(field_line_start,coord_logR,coord_theta,coord_phi,B_flicks,1.0,2.9,nlblks,n1pm1,n2pm1,n3pm1,step_size=1E-2)
		field_line_X=field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2])
		field_line_Y=-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2])
		field_line_Z=field_line_sph[:,0]*np.cos(field_line_sph[:,1])
		mlab.plot3d(field_line_X,field_line_Y,field_line_Z,line_width=0.01,color=(0,0,0),tube_radius=0.004)
		plt.plot(field_line_sph[:,2],field_line_sph[:,1],color="black")

mlab.view(azimuth=0, elevation=110, roll=90, distance=4.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)
mlab.show()



