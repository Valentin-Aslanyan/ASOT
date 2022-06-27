
"""
Manually select locations for surface flows
Double left click location to add surface flow; double right click to remove
Outputs locations in a format suitable for ARMS control files

Uses squashing factor map from QSL Squasher (not included in this distribution)
"""

phi_limits=[-30,40]	#None for default
theta_limits=[0,35]	#degrees
min_size=0.2		#degrees

v_phi=0.4E5
v_theta=-0.4E5
k_phi=40.0
k_theta=40.0
flow_width=4.0	#degrees

OCB_file="./PFLS/0000000/qslR1.bin"

#degrees
def create_limits(theta_coord,phi_coord,width):
	theta_central=np.round((90.0-point[1])/180.0*1000.0)/1000.0
	phi_central=np.round(phi_coord/180.0*1000.0)/1000.0
	width_round=np.round(width/180.0*1000)/1000.0

	theta_left=theta_central-width_round
	theta_right=theta_central+width_round
	phi_left=phi_central-width_round
	phi_right=phi_central+width_round
	return theta_left,theta_central,theta_right,phi_left,phi_central,phi_right
	


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


points_coordinates=[]
def Click_Loop(event):
	global points_coordinates
	if event.dblclick==True:
		if event.button==1 and event.xdata!=None and event.ydata!=None:	#Left click
			new_point_click=[float(event.xdata),float(event.ydata)]
			near_other_point=False
			for idx in range(len(points_coordinates)):
				if abs(points_coordinates[idx][0]-new_point_click[0])<=min_size and abs(points_coordinates[idx][1]-new_point_click[1])<=min_size:
					near_other_point=True
					break
			if near_other_point==False:
				if len(points_coordinates)==0:
					points_coordinates=[new_point_click]
				else:
					points_coordinates.append(new_point_click)

		if event.button==3 and event.xdata!=None and event.ydata!=None:	#Right click
			new_point_click=[float(event.xdata),float(event.ydata)]
			for idx in range(len(points_coordinates)-1,-1,-1):
				if abs(points_coordinates[idx][0]-new_point_click[0])<=min_size and abs(points_coordinates[idx][1]-new_point_click[1])<=min_size:
					points_coordinates.pop(idx)

	update_plot_points()
	plt.draw()


R_QSL,theta_grid_QSL,phi_grid_QSL,Q=parse_QSL_Rbinfile(OCB_file)
theta_grid_QSL=theta_grid_QSL*RAD2DEG-90.0
phi_grid_QSL*=RAD2DEG
Q_grid=np.sign(Q)*np.log(abs(Q))/np.log(10.0)
Q_grid[np.isinf(Q_grid)]=np.nan


def update_plot_points():
	if len(points_coordinates)==1:
		points_plotted.set_data([points_coordinates[0][0]],[points_coordinates[0][1]])
	elif len(points_coordinates)>2:
		points_arr=np.array(points_coordinates)
		points_plotted.set_data(points_arr[:,0],points_arr[:,1])

fig2=plt.figure("Q 2D",figsize=(8,6))
ax2=fig2.gca()
plt.title(r"$R=1R_{\odot}$",fontsize=20)
color_plot2=plt.pcolormesh(phi_grid_QSL,theta_grid_QSL,Q_grid[:-1,:-1],cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
plt.contour(phi_grid_QSL,theta_grid_QSL,np.sign(Q_grid), [0.0],colors=["black"])
points_plotted,=ax2.plot([],[],'o',color="blue")
cbar2=fig2.colorbar(color_plot2)#,ticks=[-4,-3,-2,-1,0])
cbar2.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar2.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)
plt.xlim(phi_limits)
plt.ylim(theta_limits)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)


cid_up = fig2.canvas.mpl_connect('button_press_event', Click_Loop)


plt.show()

print("      tlsfl       trsfl       tcsfl    ktsfl    vtsfl          plsfl          prsfl       pcsfl   kpsfl      vpsfl      tilsfl       tirsfl      ticsfl       ktisfl")
for point in points_coordinates:
	theta_left,theta_cnt,theta_right,phi_left,phi_cnt,phi_right=create_limits(point[1],point[0],flow_width)
	print("[      {:.3f},".format(theta_left)+"      {:.3f},".format(theta_right)+"      {:.3f},".format(theta_cnt)+"  {:+.2e},".format(k_theta)+"  {:+.2e},".format(v_theta)
+"   {:.3f},".format(phi_left)+"          {:.3f},".format(phi_right)+"       {:.3f},".format(phi_cnt)+"  {:+.2e},".format(k_phi)+"  {:+.2e},".format(v_phi)
+"    0.0e+3,       4.0e+3,    0.0e+3,     2.0],")



