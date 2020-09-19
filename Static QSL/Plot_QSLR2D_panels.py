

filenames=[
"./0050887/qslR1.bin",
"./0053109/qslR1.bin",
"./0055411/qslR1.bin",
"./0057883/qslR1.bin",
"./0060557/qslR1.bin"]

labeltext=["$t=1/16T$","$t=2/16T$","$t=3/16T$","$t=4/16T$","$t=5/16T$","$t=6/16T$","$t=7/16T$","$t=8/16T$"]

phi_limits=[-12,-2]	#None for default
theta_limits=[110,120]	#
horizontal_panels=4

plot_cage=True
cage_from_file=True
SurfaceLine_files=[
"./SurfaceLineData/SLine_4600.bin",
"./SurfaceLineData/SLine_4800.bin",
"./SurfaceLineData/SLine_5000.bin",
"./SurfaceLineData/SLine_5200.bin",
"./SurfaceLineData/SLine_5400.bin"]
CNT_file="./arms.cnt"
SurfaceLine_initial="./SurfaceLineData/SLine_4000.bin"
t_start=4000.0
delta_t=10.0
solar_Radius=1.0*7E10
plot_t=[4600.0,4800.0,5000.0,5200.0,5400.0,5600.0,5800.0,6000.0]


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


if plot_cage:
	if cage_from_file:
		static_cages,dynamic_cages_temp=surfaceline_read_binary_cage(SurfaceLine_files[0])
		dynamic_cages=[]
		for idx in range(len(dynamic_cages_temp)):
			dynamic_cages.append(np.zeros((np.shape(dynamic_cages_temp[idx])[0],np.shape(dynamic_cages_temp[idx])[1],len(SurfaceLine_files))))
			dynamic_cages[idx][:,:,0]=dynamic_cages_temp[idx][:,:]
		for idx in range(1,len(SurfaceLine_files)):
			static_cages,dynamic_cages_temp=surfaceline_read_binary_cage(SurfaceLine_files[idx])
			for idx2 in range(len(dynamic_cages_temp)):
				dynamic_cages[idx2][:,:,idx]=dynamic_cages_temp[idx2][:,:]
	else:
		static_cages,dynamic_cages_temp=surfaceline_read_binary_cage(SurfaceLine_initial)
		dynamic_cages=[]
		for idx in range(len(dynamic_cages_temp)):
			dynamic_cages.append(surfaceline_advance_cage(dynamic_cages_temp[idx],t_start,plot_t,delta_t,CNT_file,solar_Radius))


Q_list=[]
for fname in filenames:
	R,theta,phi,Q=parse_QSL_Rbinfile(fname)
	Q_grid=np.sign(Q)*np.log(abs(Q))
	Q_grid[np.isinf(Q_grid)]=0.0#np.nan		#Note: animate bug if nan in data
	Q_list.append(Q_grid)


fig=plt.figure(figsize=(15,8))
for idx in range(len(filenames)):
	plt.subplot(int(np.ceil(len(filenames)/horizontal_panels)),horizontal_panels,idx+1)
	ax=fig.gca()
	plt.tight_layout()
	color_plot=plt.pcolormesh(phi*180.0/np.pi,theta*180.0/np.pi,Q_list[idx],cmap='RdBu_r',vmin=-10,vmax=10,rasterized=True)
	if plot_cage:
		idx3=min(idx,len(plot_t)-1)
		for idx_c in range(len(dynamic_cages)):
			plt.plot(dynamic_cages[idx_c][0,:,idx3],dynamic_cages[idx_c][1,:,idx3],'-',color="black")
		for idx_c in range(len(static_cages)):
			plt.plot(static_cages[idx_c][0,:],static_cages[idx_c][1,:],color="black")

	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	if phi_limits!=None:
		plt.xlim(phi_limits)
	if theta_limits!=None:
		plt.ylim(theta_limits)
	if len(labeltext)>idx:
		if phi_limits!=None:
			if theta_limits!=None:
				plt.text(phi_limits[0]+0.1*(phi_limits[1]-phi_limits[0]),theta_limits[0]+0.8*(theta_limits[1]-theta_limits[0]),labeltext[idx],fontsize=20)
			else:
				plt.text(phi_limits[0]+0.1*(phi_limits[1]-phi_limits[0]),160,labeltext[idx],fontsize=20)
		elif theta_limits!=None:
			plt.text(-170,theta_limits[0]+0.8*(theta_limits[1]-theta_limits[0]),labeltext[idx],fontsize=20)			
		else:
			plt.text(-170,160,labeltext[idx],fontsize=20)
	if idx % horizontal_panels==0:
		plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
	else:
		ax.set_yticklabels(["","","","","","","","","","","","","",""])
	if idx>=len(filenames)-horizontal_panels:
		plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
	else:
		ax.set_xticklabels(["","","","","","","","","","","","","",""])

plt.savefig("Q_Panels.png", format="png", dpi=100,bbox_inches='tight',pad_inches=0.1)
plt.show()






