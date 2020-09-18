

solar_Radius=1.0*7E10
border_px=10
filter_px=5
filter_thresh=5
threshold_delta=0.1
reconnection_num=100
min_final_separation=1.5
max_reconnection_radius=3.0


QSL_directories=[
"./PFLS/0053109/",
"./PFLS/0054245/"]

import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
from mayavi import mlab
from scipy.interpolate import interp2d


if filter_px % 2==0:
	filter_px+=1

R1,theta_grid1,phi_grid1,Q1=parse_QSL_Rbinfile(QSL_directories[0]+'qslR1.bin')
Q_grid1=np.sign(Q1)*np.log(abs(Q1))
Q_grid1[np.isinf(Q_grid1)]=np.nan
X_grid1=np.sin(theta_grid1)*np.cos(phi_grid1)
Y_grid1=np.sin(theta_grid1)*np.sin(phi_grid1)
Z_grid1=np.cos(theta_grid1)
Q1sign_list=[np.sign(Q1)]
theta_grid1_filtered=theta_grid1[border_px:-border_px,border_px:-border_px]
phi_grid1_filtered=phi_grid1[border_px:-border_px,border_px:-border_px]
phi,theta,R,B=parse_QSL_folder(QSL_directories[0])
phi=(phi-180.0)*np.pi/180.0
theta=(theta+90.0)*np.pi/180.0
B_list=[B]
Br_list=[B[0,:,:,0]]
B_R_func0=interp2d(phi,theta,Br_list[0])

R1,theta_grid1,phi_grid1,Q1=parse_QSL_Rbinfile(QSL_directories[1]+'qslR1.bin')
Q1sign_list.append(np.sign(Q1))
reconnected=0.5-0.5*Q1sign_list[0]*Q1sign_list[1]
reconnected_filtered=np.copy(reconnected)
for idx_t in range(filter_px//2,len(reconnected[0,:])-filter_px):
	for idx_p in range(filter_px//2,len(reconnected[:,0])-filter_px):
		if sum(reconnected[idx_p-filter_px//2:idx_p+filter_px//2+1,idx_t-filter_px//2:idx_t+filter_px//2+1].flatten())<filter_thresh:
			reconnected_filtered[idx_p,idx_t]=0.0
reconnected_filtered=reconnected_filtered[border_px:-border_px,border_px:-border_px]
phi,theta,R,B=parse_QSL_folder(QSL_directories[1])
phi=(phi-180.0)*np.pi/180.0
theta=(theta+90.0)*np.pi/180.0
B_list.append(B)
Br_list.append(B[0,:,:,0])
B_R_func1=interp2d(phi,theta,Br_list[1])




connectivity=np.zeros(np.shape(reconnected_filtered))
for idx_t in range(filter_px//2,len(reconnected_filtered[0,:])-filter_px):
	for idx_p in range(filter_px//2,len(reconnected_filtered[:,0])-filter_px):
		connectivity[idx_p,idx_t]=sum(reconnected_filtered[idx_p-filter_px//2:idx_p+filter_px//2+1,idx_t-filter_px//2:idx_t+filter_px//2+1].flatten())

indices,values=find_N_max_idx(connectivity,reconnection_num)
start_indices=tidy_max_idx(indices,values,np.shape(phi_grid1_filtered)[1])


mlab.figure(bgcolor=(0.95,0.95,0.95))

sphere_mesh=mlab.mesh(X_grid1[::-1,:],Y_grid1[::-1,:],Z_grid1[:,:],scalars=Q_grid1[:,:],colormap='RdBu',vmin=-10,vmax=10)
sphere_mesh.actor.property.backface_culling = True
sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True
#mlab.colorbar(orientation="vertical")

num_lines=0
for idx in range(len(start_indices[:,0])):
	R_start=1.0
	theta_start=theta_grid1_filtered[start_indices[idx,0],start_indices[idx,1]]
	phi_start=phi_grid1_filtered[start_indices[idx,0],start_indices[idx,1]]

	field_line_start=np.array([R_start,theta_start,phi_start])
	field_line_sph1=field_line_spherical(field_line_start,R,theta,phi,B_list[0],1.0,2.9,step_size=1E-2)
	field_line_X1=field_line_sph1[:,0]*np.sin(field_line_sph1[:,1])*np.cos(field_line_sph1[:,2])
	field_line_Y1=-field_line_sph1[:,0]*np.sin(field_line_sph1[:,1])*np.sin(field_line_sph1[:,2])
	field_line_Z1=field_line_sph1[:,0]*np.cos(field_line_sph1[:,1])


	field_line_sph2=field_line_spherical(field_line_start,R,theta,phi,B_list[1],1.0,2.9,step_size=1E-2)
	field_line_X2=field_line_sph2[:,0]*np.sin(field_line_sph2[:,1])*np.cos(field_line_sph2[:,2])
	field_line_Y2=-field_line_sph2[:,0]*np.sin(field_line_sph2[:,1])*np.sin(field_line_sph2[:,2])
	field_line_Z2=field_line_sph2[:,0]*np.cos(field_line_sph2[:,1])
	
	if abs(field_line_sph2[-1,0]-field_line_sph1[-1,0])>min_final_separation and min(max(field_line_sph1[:,0]),max(field_line_sph2[:,0]))<max_reconnection_radius:
		if num_lines<1:
			mlab.plot3d(field_line_X1,field_line_Y1,field_line_Z1,line_width=0.01,color=(0,0,0),tube_radius=0.004)
			mlab.plot3d(field_line_X2,field_line_Y2,field_line_Z2,line_width=0.01,color=(0,1,0),tube_radius=0.004)
			num_lines+=1


		for idx_d in range(min(len(field_line_X1),len(field_line_X2))):
			if (field_line_X1[idx_d]-field_line_X2[idx_d])**2+(field_line_Y1[idx_d]-field_line_Y2[idx_d])**2+(field_line_Z1[idx_d]-field_line_Z2[idx_d])**2 > threshold_delta*threshold_delta:
				break
		reconnect_X=0.5*(field_line_X1[idx_d]+field_line_X2[idx_d])
		reconnect_Y=0.5*(field_line_Y1[idx_d]+field_line_Y2[idx_d])
		reconnect_Z=0.5*(field_line_Z1[idx_d]+field_line_Z2[idx_d])

		reconnect_point=mlab.points3d([reconnect_X],[reconnect_Y],[reconnect_Z],scale_factor=0.025,color=(1.0, 0.3, 0.75))
		reconnect_point.actor.property.backface_culling = True

mlab.view(azimuth=0, elevation=110, roll=90, distance=4.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)

mlab.show()

