

target_R=3.0	#
file_directory="./"
nulls_directory="./PFLS/0122766/"
flicks_file="flicks.0912893"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_file)

print('Simulation time: ',time)
print(ntblks,nlblks,newgrd)
print()
print("         Min  |  Max    grid points")
print("r:     "+"{:+.3f}".format(min(coord_logR.flatten()))+" | {:+.3f}".format(max(coord_logR.flatten())))
print("theta: "+"{:+.3f}".format(min(coord_theta.flatten()))+" | {:+.3f}".format(max(coord_theta.flatten())))
print("phi:   "+"{:+.3f}".format(min(coord_phi.flatten()))+" | {:+.3f}".format(max(coord_phi.flatten())))
print()

new_theta,new_phi,new_data=R_slice(target_R,coord_logR,coord_theta,coord_phi,data)

plt.figure("ARMS grid and field lines",figsize=(20,10))
for idx in range(len(coord_logR)):
	if np.exp(coord_logR[idx][0])<=target_R and np.exp(coord_logR[idx][1])>=target_R:
		th1=coord_theta[idx][0]/np.pi*180.0
		ph1=coord_phi[idx][0]/np.pi*180.0
		th2=coord_theta[idx][1]/np.pi*180.0
		ph2=coord_phi[idx][1]/np.pi*180.0
		plt.plot([ph1,ph1,ph2,ph2,ph1],[th1,th2,th2,th1,th1],color="grey",linewidth=1)
plt.axis('off')
plt.savefig("Grid_PhiTheta.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.figure("pflicks grid points",figsize=(8,4))
plt.title("Grid points",fontsize=20)
plt.plot(new_phi[:]/np.pi*180.0,new_theta[:]/np.pi*180.0,'o',color="red",ms=1)
plt.axis('off')
plt.xlim([-180,180])
#plt.ylim([-3,3])
plt.savefig("GridPoints_PhiTheta.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

new_rho=new_data[:,0]
new_B=np.sqrt(new_data[:,6]**2+new_data[:,7]**2+new_data[:,8]**2)
new_v=np.sqrt(new_data[:,1]**2+new_data[:,2]**2+new_data[:,3]**2)
new_J=np.sqrt(new_data[:,9]**2+new_data[:,10]**2+new_data[:,11]**2)
new_vtheta=new_data[:,2]
new_vA=new_data[:,13]

min_new_theta=min(new_theta)
max_new_theta=max(new_theta)
min_new_phi=min(new_phi)
max_new_phi=max(new_phi)

from scipy.interpolate import griddata
grid_theta=np.linspace(min_new_theta,max_new_theta,num=1600)
grid_phi=np.linspace(min_new_phi,max_new_phi,num=3200)
grid_theta,grid_phi=np.meshgrid(grid_theta,grid_phi)
grid_theta2=np.linspace(100/180.0*np.pi,120/180.0*np.pi,num=500)
grid_phi2=np.linspace(-20/180.0*np.pi,0/180.0*np.pi,num=400)
grid_theta2,grid_phi2=np.meshgrid(grid_theta2,grid_phi2)
v_grid=griddata((new_theta,new_phi), new_v, (grid_theta2,grid_phi2), method='linear')
B_grid=griddata((new_theta,new_phi), new_B, (grid_theta,grid_phi), method='linear')
J_grid=griddata((new_theta,new_phi), new_J, (grid_theta,grid_phi), method='linear')
vA_grid=griddata((new_theta,new_phi), new_vA, (grid_theta,grid_phi), method='linear')

fig1=plt.figure("pcolormesh B",figsize=(10,8))
ax1=fig1.gca()
plt.title("$B$-field,$\quad$ $R="+"{:.1f}".format(target_R)+"$",fontsize=20)
color_plot1=plt.pcolormesh(grid_phi/np.pi*180.0,grid_theta/np.pi*180.0,B_grid,cmap='bone_r',rasterized=True)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)

plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar1=fig1.colorbar(color_plot1)#,ticks=[-4,-3,-2,-1,0])
cbar1.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar1.set_label(label="$B\;[\mathrm{a.u.}]$",fontsize=20)
#plt.axis('off')
#plt.xlim([0,3.1])
#plt.ylim([-3,3])
plt.savefig("B.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

fig2=plt.figure("pcolormesh v",figsize=(5,4))
ax2=fig2.gca()
plt.title(r"$\log(|v|),\quad$ $R="+"{:.1f}".format(target_R)+"$",fontsize=20)
color_plot2=plt.pcolormesh(grid_phi2/np.pi*180.0,grid_theta2/np.pi*180.0,v_grid,cmap='Reds',rasterized=True)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar2=fig2.colorbar(color_plot2)#,ticks=[-4,-3,-2,-1,0])
cbar2.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar2.set_label(label=r"$|v|$ [km s$^{-1}$]",fontsize=20)
#plt.axis('off')
#plt.xlim([0,3.1])
#plt.ylim([-3,3])
plt.savefig("V.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

fig3=plt.figure("pcolormesh vA",figsize=(10,8))
ax3=fig3.gca()
plt.title(r"$v_A,\quad$ $R="+"{:.1f}".format(target_R)+"$",fontsize=20)
color_plot3=plt.pcolormesh(grid_phi/np.pi*180.0,grid_theta/np.pi*180.0,vA_grid,cmap='hot_r',rasterized=True)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar3=fig3.colorbar(color_plot3)#,ticks=[-4,-3,-2,-1,0])
cbar3.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar3.set_label(label=r"$v_{A}\;[\mathrm{a.u.}]$",fontsize=20)
#plt.axis('off')
#plt.xlim([0,3.1])
#plt.ylim([-3,3])
plt.savefig("V_Alfven.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


fig4=plt.figure("|J|/|B|",figsize=(14,8))
ax4=fig4.gca()
plt.title(r"$\log \left(|\vec{J}|/|\vec{B}|\right),\quad$ $R="+"{:.1f}".format(target_R)+"$",fontsize=20)
color_plot4=plt.pcolormesh(grid_phi/np.pi*180.0,grid_theta/np.pi*180.0,np.log(J_grid/B_grid)/np.log(10.0),cmap='hot_r',vmin=-4.0,vmax=0.0,rasterized=True)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar4=fig4.colorbar(color_plot4)#,ticks=[-4,-3,-2,-1,0])
cbar4.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
#cbar4.set_label(label=r"$v_{A}\;[\mathrm{a.u.}]$",fontsize=20)
#plt.axis('off')
#plt.xlim([0,3.1])
#plt.ylim([-3,3])
plt.savefig("J_B.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


plt.show()



