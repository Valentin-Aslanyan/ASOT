

t_start=0.0
t_end=8000.0

num_frames_QSL=200
num_frames_v=240
frames_per_sec=20

phi_limits=[-30,40]	#None for default
theta_limits=[90,125]	#degrees

CNT_file="./SHORT_arms.cnt"
OCB_file="./PFLS/0044295/qslR1.bin"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call


R_QSL,theta_grid_QSL,phi_grid_QSL,Q=parse_QSL_Rbinfile(OCB_file)
theta_grid_QSL*=180.0/np.pi
phi_grid_QSL*=180.0/np.pi
Q_grid=np.sign(Q)*np.log(abs(Q))/np.log(10.0)
Q_grid[np.isinf(Q_grid)]=np.nan

theta_v=np.linspace(theta_limits[0],theta_limits[1],num=300)
phi_v=np.linspace(phi_limits[0],phi_limits[1],num=400)
theta_grid_v,phi_grid_v=np.meshgrid(theta_v,phi_v)
v_grid=np.zeros((np.shape(phi_grid_v)[0],np.shape(phi_grid_v)[1],num_frames_v))


Time_profiles,sflow_parameters=parse_sflow_CNT(CNT_file)


def Curve_vortex(theta_grid,phi_grid,tlsfl,trsfl,tcsfl,ktsfl,vtsfl,plsfl,prsfl,pcsfl,kpsfl,vpsfl):
	plsfl*=np.pi
	prsfl*=np.pi
	pcsfl*=np.pi
	vpsfl*=1.0E-5
	tlsfl*=np.pi
	trsfl*=np.pi
	tcsfl*=np.pi
	vtsfl*=1.0E-5
	
	rd2rl=1.0/(trsfl-tlsfl)
	k2drl=rd2rl*ktsfl
	farg2=k2drl*rd2rl
	r_theta=1.0/np.sqrt(2.0*farg2)
	rd3rl=1.0/(prsfl-plsfl)
	k3drl=rd3rl*kpsfl
	farg3=k3drl*rd3rl
	r_phi=1.0/np.sqrt(2.0*farg3)

	chi=np.linspace(0.0,2.0*np.pi,num=200)
	theta_curve=np.pi-(r_theta*np.sin(chi)+tcsfl)
	phi_curve=r_phi*np.cos(chi)+pcsfl
	if vpsfl>0.0:
		direction=1.0
	else:
		direction=-1.0

	return phi_curve,theta_curve,direction





fig=plt.figure("v",figsize=(12,6))
ax=fig.gca()
#plt.tight_layout()
plt.title(r"$R=1R_{\odot}$",fontsize=20)
plt.contour(phi_grid_QSL,theta_grid_QSL,np.sign(Q_grid), [0.0],colors=["black"])
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
		phi_curve,theta_curve,direction=Curve_vortex(theta_grid_v,phi_grid_v,tlsfl,trsfl,tcsfl,ktsfl,vtsfl,plsfl,prsfl,pcsfl,kpsfl,vpsfl)
		if direction>0.0:
			plot_color="orange"
		else:
			plot_color="green"
		plt.plot(phi_curve*180.0/np.pi,theta_curve*180.0/np.pi,color=plot_color,linewidth=3)

#cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
#cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
#cbar.set_label(label=r"$|v|$ [km s$^{-1}$]",fontsize=20)
plt.xlim(phi_limits)
plt.ylim(theta_limits)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)

plt.savefig("SurfaceVelocityDirection.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.show()



