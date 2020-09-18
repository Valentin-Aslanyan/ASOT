

R=1.0*7E10
phi_lim=(-25,10)
theta_lim=(95,130)
cage_dimension=(10,5)
cage_phi_lim=(-20,5)
cage_theta_lim=(114,125)
ignore_uniform_velocity=True
CNT_file="./arms.cnt"
t_start=4.295e+3
t_end=1.35e+4
delta_t=10.0
plot_t=[4300.0,4800.0,5300.0,5700.0,6000.0,6300.0,6700.0,7100.0,7500.0]

frames_per_step=3
frames_per_sec=3
pad_start_frames=2
pad_end_frames=2


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call


Time_profiles,sflow_parameters=parse_sflow_CNT(CNT_file)

######Set up cage
t=t_start
num_t=int((t_end-t_start)/delta_t+1)
cage=np.zeros((2,cage_dimension[0]*cage_dimension[1]))
cage_output=np.zeros((2,cage_dimension[0]*cage_dimension[1],len(plot_t)))
cage_connections=np.zeros((2,cage_dimension[0]*(cage_dimension[1]-1)+(cage_dimension[0]-1)*cage_dimension[1]),dtype=np.int32)
idx_cc=0
for idx_t in range(cage_dimension[1]):
	for idx_p in range(cage_dimension[0]):
		cdim1m_0=cage_dimension[0]
		if cdim1m_0==0:
			cdim1m_0=1
		cdim1m_1=cage_dimension[1]
		if cdim1m_1==0:
			cdim1m_1=1
		cage[0,idx_t*cage_dimension[0]+idx_p]=(cage_phi_lim[0]+(cage_phi_lim[1]-cage_phi_lim[0])*idx_p/cdim1m_0)/180.0
		cage[1,idx_t*cage_dimension[0]+idx_p]=1.0-(cage_theta_lim[0]+(cage_theta_lim[1]-cage_theta_lim[0])*idx_t/cdim1m_1)/180.0
		if idx_t>0:
			cage_connections[0,idx_cc]=idx_t*cage_dimension[0]+idx_p
			cage_connections[1,idx_cc]=(idx_t-1)*cage_dimension[0]+idx_p
			idx_cc+=1
		if idx_p>0:
			cage_connections[0,idx_cc]=idx_t*cage_dimension[0]+idx_p
			cage_connections[1,idx_cc]=idx_t*cage_dimension[0]+(idx_p-1)
			idx_cc+=1
#########################

######Evolve cage
idx_t_out=0
for idx_t in range(1,num_t):
	for idx_c in range(len(cage[0,:])):
		new_pos=surface_move_rk_4(R,cage[1,idx_c],cage[0,idx_c],t,delta_t,Time_profiles,sflow_parameters)
		cage[0,idx_c]=new_pos[1]
		cage[1,idx_c]=new_pos[0]
	if idx_t_out<len(plot_t) and t>plot_t[idx_t_out]:
		cage_output[:,:,idx_t_out]=cage[:,:]
		idx_t_out+=1
	t+=delta_t
#################

call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(plot_t)):
	plt.clf()
	fig=plt.figure(figsize=(8,8))
	ax=fig.gca()
	plt.tight_layout()
	plt.xlim(phi_lim)
	plt.ylim(theta_lim)
	plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
	plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	plt.plot(cage_output[0,:,idx_f]*180.0,(1.0-cage_output[1,:,idx_f])*180.0,'o',color="blue")

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
	if idx_f==len(plot_t)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_SurfacePoints.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



