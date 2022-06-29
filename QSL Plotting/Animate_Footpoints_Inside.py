
"""
Display squashing factor from QSL Squasher (not included in this distribution) at, and above, the photosphere
Display a surface line or "cage" (corresponding to e.g. open closed boundary or Polarity Inversion Line) on the photosphere plot
Draw a locus of the surface line at multiple given distances away and identify footpoints of fieldlines along it; trace the fieldlines upwards to "handholds" in the above plot
SurfaceLine_file must have single static cage

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""

phi_points=400
phi_limits=[-30.0,40.0]
phi_limits_outer=[-65,65]
theta_points=200
theta_limits=[0.0,35.0]
theta_limits_outer=[-20,50]
num_footpoints=20
file_directory="./PFLS/0000000/"
contour_distances=[0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0]
SurfaceLine_file="./SurfaceLineData/SLine_StaticStart.bin"

frames_per_step=3
frames_per_sec=2
pad_start_frames=3
pad_end_frames=3


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call
from scipy.interpolate import interp1d


R1,theta_grid1,phi_grid1,Q1=parse_QSL_Rbinfile(file_directory+'qslR1.bin')
Q_grid1=np.sign(Q1)*np.log(abs(Q1))
Q_grid1[np.isinf(Q_grid1)]=np.nan


R3,theta_grid3,phi_grid3,Q3=parse_QSL_Rbinfile(file_directory+'qslR3.bin')
Q_grid3=np.sign(Q3)*np.log(abs(Q3))
Q_grid3[np.isinf(Q_grid3)]=np.nan

phi,theta,R,B=parse_QSL_folder(file_directory)
phi=(phi-180.0)*np.pi/180.0
theta=(theta+90.0)*np.pi/180.0

static_cages,dynamic_cages=surfaceline_read_binary_cage(SurfaceLine_file)

phi_grid=np.linspace(phi_limits[0],phi_limits[1],num=phi_points)
theta_grid=np.linspace(theta_limits[0]+90.0,theta_limits[1]+90.0,num=theta_points)

distance2=np.zeros((theta_points,phi_points))
for idx_p in range(phi_points):
	for idx_t in range(theta_points):
		cage_distance2=(static_cages[0][0,:]-phi_grid[idx_p])**2+(static_cages[0][1,:]-theta_grid[idx_t])**2
		distance2[idx_t,idx_p]=min(cage_distance2)

phi_grid,theta_grid=np.meshgrid(phi_grid,theta_grid)

footpoints=np.zeros((len(contour_distances),2,num_footpoints))
handholds=np.zeros((len(contour_distances),2,num_footpoints))
min_th=min(static_cages[0][1,:])
max_th=max(static_cages[0][1,:])
min_ph=min(static_cages[0][0,:])
max_ph=max(static_cages[0][0,:])
inner_contour=[]
for idx_f in range(len(contour_distances)):
	fig=plt.figure()
	contour1=plt.contour(phi_grid,theta_grid,distance2,[contour_distances[idx_f]**2],colors=["grey"])
	plt.close(fig)

	for idx in range(len(contour1.allsegs[0])):
		points_inside=True
		for idx2 in range(len(contour1.allsegs[0][idx][:,0])):
			if contour1.allsegs[0][idx][idx2,0]>max_ph or contour1.allsegs[0][idx][idx2,0]<min_ph or contour1.allsegs[0][idx][idx2,1]>max_th or contour1.allsegs[0][idx][idx2,1]<min_th:
				points_inside=False
				break
		if points_inside:
			inner_contour.append(np.transpose(contour1.allsegs[0][idx]))
			break

	x_param=np.linspace(1,num_footpoints,num=len(inner_contour[idx_f][0,:]))
	x_out=np.linspace(1,num_footpoints,num=num_footpoints)
	phi_func=interp1d(x_param,inner_contour[idx_f][0,:],kind='linear')
	theta_func=interp1d(x_param,inner_contour[idx_f][1,:],kind='linear')

	footpoints[idx_f,0,:]=phi_func(x_out)
	footpoints[idx_f,1,:]=theta_func(x_out)


	for idx in range(num_footpoints):
		field_line_start=np.array([1.0,footpoints[idx_f,1,idx]*DEG2RAD,footpoints[idx_f,0,idx]*DEG2RAD])
		field_line_sph=field_line_spherical(field_line_start,R,theta,phi,B,1.0,3.0,step_size=1E-2)
		handholds[idx_f,0,idx]=field_line_sph[-1,2]*RAD2DEG
		handholds[idx_f,1,idx]=field_line_sph[-1,1]*RAD2DEG

call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(contour_distances)):
	plt.clf()
	fig=plt.figure("Qmap and points",figsize=(12,13))
	plt.subplot(2,1,1)
	plt.text(phi_limits[0]+2,theta_limits[0]+2,r"$R=1R_{\odot}$",fontsize=26)
	colour_plot=plt.pcolormesh(phi_grid1*RAD2DEG,theta_grid1*RAD2DEG-90.0,Q_grid1[:-1,:-1],cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
	plt.plot(inner_contour[idx_f][0,:],inner_contour[idx_f][1,:]-90.0,color="grey")
	plt.plot(static_cages[0][0,:],static_cages[0][1,:]-90.0,color="black")
	plt.plot(footpoints[idx_f,0,:],footpoints[idx_f,1,:]-90.0,'o',color="green")

	plt.xlim(phi_limits)
	plt.ylim(theta_limits)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
	plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)

	plt.subplot(2,1,2)
	plt.text(phi_limits_outer[0]+3,theta_limits_outer[1]-11,r"$R=3R_{\odot}$",fontsize=26)
	colour_plot=plt.pcolormesh(phi_grid3*RAD2DEG,theta_grid3*RAD2DEG-90.0,Q_grid3[:-1,:-1],cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
	plt.plot(handholds[idx_f,0,:],handholds[idx_f,1,:]-90.0,'o',color="green")

	plt.xlim(phi_limits_outer)
	plt.ylim(theta_limits_outer)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
	plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)


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
	if idx_f==len(contour_distances)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Footpoints_Inside.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



