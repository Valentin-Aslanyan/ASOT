
"""
Animate "cages" (e.g. open/closed boundary or Polarity Inversion Line) on top of squashing factor from QSL Squasher (not included in this distribution)
"""

solar_Radius=1.0*7E10
phi_limits=[-30,40]	#None for default
theta_limits=[0,35]	#degrees

frames_per_step=3
frames_per_sec=2
pad_start_frames=6
pad_end_frames=6


cage_filenames=[
"./SurfaceLineData/SLine_1000.bin",
"./SurfaceLineData/SLine_2000.bin",
"./SurfaceLineData/SLine_3000.bin",
"./SurfaceLineData/SLine_4000.bin",
"./SurfaceLineData/SLine_5000.bin"]

QSL_filenames=[
"./PFLS/0001000/qslR1.bin",
"./PFLS/0002000/qslR1.bin",
"./PFLS/0003000/qslR1.bin",
"./PFLS/0004000/qslR1.bin",
"./PFLS/0005000/qslR1.bin"]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from subprocess import call


Q_list=[]
for fname in QSL_filenames:
	R,theta,phi,Q=parse_QSL_Rbinfile(fname)
	Q_grid=np.sign(Q)*np.log(abs(Q))/np.log(10.0)
	Q_grid[np.isinf(Q_grid)]=0.0#np.nan		#Note: animate bug if nan in data
	Q_list.append(Q_grid)


static_cages=[]
dynamic_cages=[]
for fname in cage_filenames:
	sc,dc=surfaceline_read_binary_cage(cage_filenames)
	static_cages.append(sc)
	dynamic_cages.append(dc)


call_result=call(["mkdir","./anim_temp"])


plot_idx=0
for idx_f in range(min(len(cage_filenames),len(QSL_filenames))):
	plt.clf()
	fig=plt.figure("out",figsize=(8,8))
	ax=fig.gca()

	plt.tight_layout()
	plt.xlim(phi_limits)
	plt.ylim(theta_limits)
	plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
	plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)

	color_plot=plt.pcolormesh(phi*RAD2DEG,theta*RAD2DEG-90.0,Q_list[idx_f],cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
	for idx in range(len(static_cages[idx_f])):
		plt.plot(static_cages[idx_f][idx][0,:],static_cages[idx_f][idx][1,:]-90.0,color="black")
	for idx in range(len(dynamic_cages[idx_f])):
		plt.plot(dynamic_cages[idx_f][idx][0,:],dynamic_cages[idx_f][idx][1,:]-90.0,color="black")

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
	if idx_f==min(len(cage_filenames),len(QSL_filenames))-1:
		for idx in range(pad_end_frames):
			call_result=call(["cp","./anim_temp/img{:03d}.png".format(plot_idx-1),"./anim_temp/img{:03d}.png".format(plot_idx)])
			plot_idx+=1


call_result=call(['ffmpeg -framerate '+str(frames_per_sec)+' -i ./anim_temp/img%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ./anim_temp/anim.mp4'],shell=True)
call_result=call(["cp","./anim_temp/anim.mp4","./Animate_SurfaceLine.mp4"])
call_result=call(["rm","-r","./anim_temp/"])



