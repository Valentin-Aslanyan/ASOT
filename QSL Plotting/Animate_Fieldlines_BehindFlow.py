
"""
3D sphere with squashing factor from QSL Squasher (not included in this distribution) at, and optionally above, the photosphere
Display a surface line or "cage" (corresponding to e.g. open closed boundary or Polarity Inversion Line) on the lower sphere
Draw specified field lines upwards from the photosphere

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""

QSL_directories=[
"./PFLS/0000001/",
"./PFLS/0000002/",
"./PFLS/0000003/",
"./PFLS/0000004/",
"./PFLS/0000005/"]

SurfaceLine_filenames=[
",/SurfaceLineData/SLine_1000.bin",
",/SurfaceLineData/SLine_2000.bin",
",/SurfaceLineData/SLine_3000.bin",
",/SurfaceLineData/SLine_4000.bin",
",/SurfaceLineData/SLine_5000.bin"]

R_start    =[ +1.000, +1.000, +1.000, +1.000, +1.000, +1.000]
theta_start=[+24.477,+25.107,+24.248,+22.586,+22.930,+23.388]
phi_start  =[ -7.391, -8.251, -9.339, -7.105, -7.907, -5.730]
plot_R3=False

frames_per_step=3
frames_per_sec=2
pad_start_frames=6
pad_end_frames=6


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
from mayavi import mlab
from subprocess import call


num_timesteps=min(len(QSL_directories),len(SurfaceLine_filenames))

Q1_list=[]
Q3_list=[]
for idx in range(num_timesteps):
	R1,theta_grid1,phi_grid1,Q1=parse_QSL_Rbinfile(QSL_directories[idx]+'qslR1.bin')
	Q_grid1=np.sign(Q1)*np.log(abs(Q1))
	Q_grid1[np.isinf(Q_grid1)]=0.0#np.nan		#Note: animate bug if nan in data
	Q1_list.append(Q_grid1)
	if plot_R3:
		R3,theta_grid3,phi_grid3,Q3=parse_QSL_Rbinfile(QSL_directories[idx]+'qslR3.bin')
		Q_grid3=np.sign(Q3)*np.log(abs(Q3))
		Q_grid3[np.isinf(Q_grid3)]=0.0
		Q3_list.append(Q_grid3)

X_grid1=np.sin(theta_grid1)*np.cos(phi_grid1)
Y_grid1=np.sin(theta_grid1)*np.sin(phi_grid1)
Z_grid1=np.cos(theta_grid1)

if plot_R3:
	X_grid3=2.89*np.sin(theta_grid3)*np.cos(phi_grid3)
	Y_grid3=2.89*np.sin(theta_grid3)*np.sin(phi_grid3)
	Z_grid3=2.89*np.cos(theta_grid3)


num_fieldlines=min(len(R_start),len(theta_start),len(phi_start))

field_line_X=[]
field_line_Y=[]
field_line_Z=[]
for idx in range(num_timesteps):
	phi,theta,R,B=parse_QSL_folder(QSL_directories[idx])
	phi=(phi-180.0)*np.pi/180.0
	theta=(theta+90.0)*np.pi/180.0

	field_line_X.append([])
	field_line_Y.append([])
	field_line_Z.append([])
	for idx2 in range(num_fieldlines):
		th_temp,ph_temp=change_angular_coords(theta_start[idx2],phi_start[idx2],from_type='carrington',to_type='flicks')
		field_line_start=np.array([R_start[idx2],th_temp+0.5*np.pi,ph_temp])
		field_line_sph=field_line_spherical(field_line_start,R,theta,phi,B,1.0,2.9,step_size=1E-2)
		field_line_X[idx].append(field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2]))
		field_line_Y[idx].append(-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2]))
		field_line_Z[idx].append(field_line_sph[:,0]*np.cos(field_line_sph[:,1]))


static_cages=[]
dynamic_cages=[]
for idx_f in range(num_timesteps):
	static_cages_temp,dynamic_cages_temp=surfaceline_read_binary_cage(SurfaceLine_filenames[idx_f])
	static_cages.append(static_cages_temp)
	dynamic_cages.append(dynamic_cages_temp)

call_result=call(["mkdir","./anim_temp"])


#Do the 0th frame twice, first one has a problem
mlab.figure(bgcolor=(0.0,0.0,0.0),size=(1200,1000))

mlab.clf()

sphere_mesh=mlab.mesh(X_grid1[::-1,:],Y_grid1[::-1,:],Z_grid1[:,:],scalars=Q1_list[0][:,:],colormap='RdBu',vmin=-10,vmax=10)
sphere_mesh.actor.property.backface_culling = True
sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True
#mlab.colorbar(orientation="vertical")

if plot_R3:
	sphere_mesh3=mlab.mesh(X_grid3[::-1,:],Y_grid3[::-1,:],Z_grid3[:,:],scalars=Q3_list[0][:,:],colormap='RdBu',vmin=-10,vmax=10,opacity=0.4)
	sphere_mesh3.actor.property.backface_culling = True
	sphere_mesh3.module_manager.scalar_lut_manager.reverse_lut = True

for idx_c in range(len(dynamic_cages[0])):
	cage_X=1.0*np.sin(dynamic_cages[0][idx_c][1,:]*DEG2RAD)*np.cos(dynamic_cages[0][idx_c][0,:]*DEG2RAD)
	cage_Y=-1.0*np.sin(dynamic_cages[0][idx_c][1,:]*DEG2RAD)*np.sin(dynamic_cages[0][idx_c][0,:]*DEG2RAD)
	cage_Z=1.0*np.cos(dynamic_cages[0][idx_c][1,:]*DEG2RAD)
	mlab.plot3d(cage_X,cage_Y,cage_Z,line_width=0.01,color=(0,0,0),tube_radius=0.001)
for idx_c in range(len(static_cages[0])):
	cage_X=1.0*np.sin(static_cages[0][idx_c][1,:]*DEG2RAD)*np.cos(static_cages[0][idx_c][0,:]*DEG2RAD)
	cage_Y=-1.0*np.sin(static_cages[0][idx_c][1,:]*DEG2RAD)*np.sin(static_cages[0][idx_c][0,:]*DEG2RAD)
	cage_Z=1.0*np.cos(static_cages[0][idx_c][1,:]*DEG2RAD)
	mlab.plot3d(cage_X,cage_Y,cage_Z,line_width=0.01,color=(0,0,0),tube_radius=0.001)


#Fudge to set figure limits
mlab.plot3d([2.9,3.0],[3.0,3.0],[3.0,3.0],line_width=0.01,color=(0,0,0),tube_radius=0.001)
mlab.plot3d([-2.9,-3.0],[-3.0,-3.0],[-3.0,-3.0],line_width=0.01,color=(0,0,0),tube_radius=0.001)

mlab.view(azimuth=-20, elevation=110, roll=90, distance=5.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)

mlab.savefig("./anim_temp/img{:03d}.png".format(0))

fig=mlab.gcf()
mlab.savefig("./anim_temp/img{:03d}.png".format(0))
fig.scene.light_manager.lights[0].azimuth=-5.0
fig.scene.light_manager.lights[0].intensity=0.8

#0th frame properly
plot_idx=0
for idx_f in range(num_timesteps):
	mlab.clf()
	fig=mlab.gcf()

	sphere_mesh=mlab.mesh(X_grid1[::-1,:],Y_grid1[::-1,:],Z_grid1[:,:],scalars=Q1_list[idx_f][:,:],colormap='RdBu',vmin=-10,vmax=10)
	sphere_mesh.actor.property.backface_culling = True
	sphere_mesh.module_manager.scalar_lut_manager.reverse_lut = True
	#mlab.colorbar(orientation="vertical")

	for idx_c in range(len(dynamic_cages[idx_f])):
		cage_X=1.0*np.sin(dynamic_cages[idx_f][idx_c][1,:]*DEG2RAD)*np.cos(dynamic_cages[idx_f][idx_c][0,:]*DEG2RAD)
		cage_Y=-1.0*np.sin(dynamic_cages[idx_f][idx_c][1,:]*DEG2RAD)*np.sin(dynamic_cages[idx_f][idx_c][0,:]*DEG2RAD)
		cage_Z=1.0*np.cos(dynamic_cages[idx_f][idx_c][1,:]*DEG2RAD)
		mlab.plot3d(cage_X,cage_Y,cage_Z,line_width=0.01,color=(0,0,0),tube_radius=0.001)
	for idx_c in range(len(static_cages[idx_f])):
		cage_X=1.0*np.sin(static_cages[idx_f][idx_c][1,:]*DEG2RAD)*np.cos(static_cages[idx_f][idx_c][0,:]*DEG2RAD)
		cage_Y=-1.0*np.sin(static_cages[idx_f][idx_c][1,:]*DEG2RAD)*np.sin(static_cages[idx_f][idx_c][0,:]*DEG2RAD)
		cage_Z=1.0*np.cos(static_cages[idx_f][idx_c][1,:]*DEG2RAD)
		mlab.plot3d(cage_X,cage_Y,cage_Z,line_width=0.01,color=(0,0,0),tube_radius=0.001)

	#Fudge to set figure limits
	mlab.plot3d([2.9,3.0],[3.0,3.0],[3.0,3.0],line_width=0.01,color=(0,0,0),tube_radius=0.001)
	mlab.plot3d([-2.9,-3.0],[-3.0,-3.0],[-3.0,-3.0],line_width=0.01,color=(0,0,0),tube_radius=0.001)


	if plot_R3:
		sphere_mesh3=mlab.mesh(X_grid3[::-1,:],Y_grid3[::-1,:],Z_grid3[:,:],scalars=Q3_list[idx_f][:,:],colormap='RdBu',vmin=-10,vmax=10,opacity=0.4)
		sphere_mesh3.actor.property.backface_culling = True
		sphere_mesh3.module_manager.scalar_lut_manager.reverse_lut = True

	min_X=min(field_line_X[idx_f][0])
	max_X=max(field_line_X[idx_f][0])
	min_Y=min(field_line_Y[idx_f][0])
	max_Y=max(field_line_Y[idx_f][0])
	min_Z=min(field_line_Z[idx_f][0])
	max_Z=max(field_line_Z[idx_f][0])
	for idx_l in range(num_fieldlines):
		mlab.plot3d(field_line_X[idx_f][idx_l],field_line_Y[idx_f][idx_l],field_line_Z[idx_f][idx_l],line_width=0.01,color=(237/255,166/255,0),tube_radius=0.001)
		min_X=min(min(field_line_X[idx_f][idx_l]),min_X)
		max_X=max(max(field_line_X[idx_f][idx_l]),max_X)
		min_Y=min(min(field_line_Y[idx_f][idx_l]),min_X)
		max_Y=max(max(field_line_Y[idx_f][idx_l]),max_Y)
		min_Z=min(min(field_line_Z[idx_f][idx_l]),min_X)
		max_Z=max(max(field_line_Z[idx_f][idx_l]),max_Z)

	mlab.savefig("./anim_temp/img{:03d}.png".format(idx_f))
	mlab.view(azimuth=-20, elevation=110, roll=90, distance=5.0)#, focalpoint=None, roll=None, reset_roll=True, figure=None)
	fig.scene.light_manager.lights[0].azimuth=-5.0
	fig.scene.light_manager.lights[0].intensity=0.8

	mlab.savefig("./anim_temp/img{:03d}.png".format(idx_f))
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(1,frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==num_timesteps-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Fieldlines_BehindFlow.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



