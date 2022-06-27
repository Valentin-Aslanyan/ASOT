
"""
Display time-dependent map of |v| for surface flow vortices (directly from ARMS control file) relative to Open/Closed boundary (or other "cage", e.g. Polarity Inversion Line)
"""

t_start=0.0
t_end=1000.0

num_frames_v=240
frames_per_step=3
frames_per_sec=10
pad_start_frames=3
pad_end_frames=3

phi_limits=[-30,40]	#None for default
theta_limits=[0,35]	#degrees

CNT_file="./arms.cnt"
cage_file="./SLine.bin"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call


size_t_grid=400
size_p_grid=500


theta_v=np.linspace(theta_limits[0]+90.0,theta_limits[1]+90.0,num=size_t_grid)
phi_v=np.linspace(phi_limits[0],phi_limits[1],num=size_p_grid)
theta_grid_v,phi_grid_v=np.meshgrid(theta_v,phi_v)
v_grid=np.zeros((np.shape(phi_grid_v)[0],np.shape(phi_grid_v)[1],num_frames_v))


Time_profiles,sflow_parameters=parse_sflow_CNT(CNT_file)
static_cages,dynamic_cages=surfaceline_read_binary_cage(cage_file)


def GaussDGauss_vortex(theta_grid,phi_grid,tlsfl,trsfl,tcsfl,ktsfl,vtsfl,plsfl,prsfl,pcsfl,kpsfl,vpsfl):
	plsfl*=np.pi
	prsfl*=np.pi
	pcsfl*=np.pi
	vpsfl*=1.0E-5
	tlsfl*=np.pi
	trsfl*=np.pi
	tcsfl*=np.pi
	vtsfl*=1.0E-5

	theta_grid_new=(180.0-theta_grid)/180.0*np.pi-tcsfl
	phi_grid_new=phi_grid/180.0*np.pi-pcsfl
	
	rd2rl=1.0/(trsfl-tlsfl)
	k2drl=rd2rl*ktsfl
	rd3rl=1.0/(prsfl-plsfl)
	k3drl=rd3rl*kpsfl
	farg3=k3drl*rd3rl*phi_grid_new
	garg3=farg3*phi_grid_new
	gexp3=np.exp(-garg3)
	farg2=k2drl*rd2rl*theta_grid_new
	garg2=farg2*theta_grid_new
	gexp2=np.exp(-garg2)

	v_theta=-2.0*vtsfl*farg3*gexp3*gexp2/np.sin(np.pi-theta_grid/180.0*np.pi)
	v_phi=-2.0*vpsfl*farg2*gexp2*gexp3

	return v_theta,v_phi

def ARMS_cos_t(t,tilsfl,tirsfl,ticsfl,ktisfl):
	if t>tirsfl or t<tilsfl:
		cos_t=0.0
	else:
		t_p=2.0*np.pi*ktisfl/(tirsfl-tilsfl)
		cos_t=0.5*(1.0-np.cos(t_p*(t-ticsfl)))
	return cos_t


for idx_t in range(num_frames_v):
	v_theta_tot=np.zeros(np.shape(theta_grid_v))
	v_phi_tot=np.zeros(np.shape(theta_grid_v))
	t_current=t_start+(t_end-t_start)*idx_t/(num_frames_v-1)
	for idx in range(len(sflow_parameters)):
		if Time_profiles[idx]>1:
			tlsfl=sflow_parameters[idx][5]
			trsfl=sflow_parameters[idx][6]
			tcsfl=sflow_parameters[idx][7]
			ktsfl=sflow_parameters[idx][8]
			vtsfl=sflow_parameters[idx][9]
			plsfl=sflow_parameters[idx][10]
			prsfl=sflow_parameters[idx][11]
			pcsfl=sflow_parameters[idx][12]
			kpsfl=sflow_parameters[idx][13]
			vpsfl=sflow_parameters[idx][14]
			tilsfl=sflow_parameters[idx][15]
			tirsfl=sflow_parameters[idx][16]
			ticsfl=sflow_parameters[idx][17]
			ktisfl=sflow_parameters[idx][18]
			v_theta,v_phi=GaussDGauss_vortex(theta_grid_v,phi_grid_v,tlsfl,trsfl,tcsfl,ktsfl,vtsfl,plsfl,prsfl,pcsfl,kpsfl,vpsfl)
			t_func=ARMS_cos_t(t_current,tilsfl,tirsfl,ticsfl,ktisfl)
			v_theta_tot+=v_theta*t_func
			v_phi_tot+=v_phi*t_func
	v_grid[:,:,idx_t]=np.sqrt(v_theta_tot**2+v_phi_tot**2)

max_v=max(v_grid.flatten())

call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(num_frames_v):
	plt.clf()
	fig=plt.figure("v",figsize=(8,6))
	ax=fig.gca()
	#plt.tight_layout()
	plt.title(r"$R=1R_{\odot}$",fontsize=20)
	color_plot=plt.pcolormesh(phi_grid_v,theta_grid_v-90.0,v_grid[:-1,:-1,idx_f],cmap='OrRd',vmin=0.0,vmax=max_v)
	for idx in range(len(static_cages)):
		plt.plot(static_cages[idx][0,:],static_cages[idx][1,:]-90.0,color="black")
	for idx in range(len(dynamic_cages)):
		plt.plot(dynamic_cages[idx][0,:],dynamic_cages[idx][1,:]-90.0,color="black")
	cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
	cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
	cbar.set_label(label=r"$|v|$ [km s$^{-1}$]",fontsize=20)
	plt.xlim(phi_limits)
	plt.ylim(theta_limits)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
	plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)

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
	if idx_f==num_frames_v-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./SurfaceVelocityMagnitude.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



