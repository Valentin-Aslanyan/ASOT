
"""
Animate successive Q maps from QSLsquasher (modified version outputting binary data not included in this repository)

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""


filenames=[
"./PFLS/0000000/qslR1.bin",
"./PFLS/0000001/qslR1.bin",
"./PFLS/0000002/qslR1.bin",
"./PFLS/0000003/qslR1.bin",
"./PFLS/0000004/qslR1.bin",
"./PFLS/0000005/qslR1.bin"]

frames_per_step=4
frames_per_sec=2
pad_start_frames=2
pad_end_frames=2

phi_limits=[-50,50]	#None for default
theta_limits=[10,60]	#


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call

Q_list=[]
for fname in filenames:
	R,theta,phi,Q=parse_QSL_Rbinfile(fname)
	Q_grid=np.sign(Q)*np.log(abs(Q))/np.log(10.0)
	Q_grid[np.isinf(Q_grid)]=0.0#np.nan		#Note: animate bug if nan in data
	Q_list.append(Q_grid)


call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(filenames)):
	plt.clf()
	fig=plt.figure(figsize=(13,6))
	ax=fig.gca()
	plt.tight_layout()
	color_plot=plt.pcolormesh(phi*RAD2DEG,theta*RAD2DEG-90.0,Q_list[idx_f],cmap='RdBu_r',vmin=-5,vmax=5)
	plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
	plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	if phi_limits!=None:
		plt.xlim(phi_limits)
	if theta_limits!=None:
		plt.ylim(theta_limits)

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
	if idx_f==len(filenames)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_QSLR2D.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



