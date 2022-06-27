
"""
Animate successive Q maps from fieldlines (not included in this repository)

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""


target_files=[
"fieldline.0000000",
"fieldline.0000001",
"fieldline.0000002",
"fieldline.0000003",
"fieldline.0000004"]

frames_per_step=3
frames_per_sec=10
pad_start_frames=3
pad_end_frames=3


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call



call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(len(target_files)):
	R,theta,phi,B_ends,flines,mapping,Q=read_fieldline_file(target_files[idx_f])
	phi_grid,theta_grid=np.meshgrid(phi,theta)
	max_Q=max(abs(Q).flatten())
	Q_grid=np.sign(Q)*np.log(np.clip(abs(Q),2.0,max_Q))/np.log(10.0)
	Q_grid[np.isinf(Q_grid)]=np.nan	

	plt.clf()
	fig=plt.figure("Output",figsize=(20,10))
	color_plot=plt.pcolormesh(phi_grid,theta_grid,Q_grid,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
	plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
	cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
	cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
	cbar.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)

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
	if idx_f==len(target_files)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Fieldlines.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



