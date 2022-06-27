
"""
Animate successive connectivity maps, see Connectivity_Map_Legend.py for colors

Requires ffmpeg; written for Linux machines (possibly MacOS), needs altering for Windows
pad_start_frames and pad_end_frames will repeat first/last frames
"""


connection_filenames=[
"./Connection_t0.bin",
"./Connection_t6000.bin",
"./Connection_t12000.bin",
"./Connection_t18000.bin",
"./Connection_t24000.bin",
"./Connection_t30000.bin",
"./Connection_t36000.bin"]

titles=[
r"$R=20R_{\odot}\quad$ $t=\;\;\;\;0$ s",
r"$R=20R_{\odot}\quad$ $t=\;6000$ s",
r"$R=20R_{\odot}\quad$ $t=12000$ s",
r"$R=20R_{\odot}\quad$ $t=18000$ s",
r"$R=20R_{\odot}\quad$ $t=24000$ s",
r"$R=20R_{\odot}\quad$ $t=30000$ s",
r"$R=20R_{\odot}\quad$ $t=36000$ s"]

phi_limits=[-70,70]
theta_limits=[-35,35]

plot_PIL=True
PIL_cagefiles=[
"./SLine_PILr20_0.bin",
"./SLine_PILr20_6000.bin",
"./SLine_PILr20_12000.bin",
"./SLine_PILr20_18000.bin",
"./SLine_PILr20_24000.bin",
"./SLine_PILr20_30000.bin",
"./SLine_PILr20_36000.bin"]

frames_per_step=2
frames_per_sec=2
pad_start_frames=2
pad_end_frames=2


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call


connection_colors=[(171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)]
connection_cmap=matplotlib.colors.ListedColormap(((171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)))


call_result=call(["mkdir","./anim_temp"])

plot_idx=0
for idx_f in range(len(connection_filenames)):

	theta_grid,phi_grid,connection_map=load_connection_map(connection_filenames[idx_f])
	connection_map=np.floor(connection_map/2)

	theta_spacing=theta_grid[0,1]-theta_grid[0,0]
	phi_spacing=phi_grid[1,0]-phi_grid[0,0]

	plt.clf()
	fig1=plt.figure("Connection map",figsize=(16,10))
	ax1_1=fig1.gca()
	color_plot1_1=plt.pcolormesh((phi_grid-0.5*phi_spacing)/np.pi*180.0,(theta_grid-0.5*theta_spacing)/np.pi*180.0-90.0,connection_map+0.5,cmap=connection_cmap,vmin=0,vmax=5,rasterized=True)
	if plot_PIL:
		static_cages,dynamic_cages=surfaceline_read_binary_cage(PIL_cagefiles[idx_f])
		for idx in range(len(static_cages)):
			plt.plot(static_cages[idx][0,:],static_cages[idx][1,:]-90.0,color="black")
	plt.title(titles[idx_f],fontsize=46)
	plt.tick_params(axis='both', which='major',labelsize=40,direction='in',bottom=True, top=True, left=True, right=True)
	plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=42)
	plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=42)
	plt.xlim(phi_limits)
	plt.ylim(theta_limits)

	plt.savefig("./anim_temp/img{:03d}.png".format(plot_idx), format="png", dpi=100,bbox_inches='tight',pad_inches=0.1)
	plt.close(fig1)
	plot_idx+=1
	if idx_f==0:
		for idx in range(pad_start_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1
	for idx in range(frames_per_step):
		call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
		plot_idx+=1
	if idx_f==len(connection_filenames)-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_Connectivity_Map.mp4"])
call_result=call(["rm","-r","./anim_temp/"])







