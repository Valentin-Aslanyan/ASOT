
"""
Display expansion factor (B_upper/B_lower, where upper/lower are taken at two radii) from modified version of QSL Squasher (not included in this distribution)
"""

filename="./PFLS/0000000/qslR3Expansion.bin"
theta_limits=[-70,60]
phi_limits=[-100,100]

plot_contour=False


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')



R_actual,grid_theta,grid_phi,Q_grid,Bs_grid,Be_grid=parse_QSL_RbinfileExpansion(filename)
Expansion_grid=np.zeros(np.shape(Q_grid))
for idx_t in range(len(grid_theta[0,:])):
	for idx_p in range(len(grid_phi[:,0])):
		if not np.isinf(Q_grid[idx_p,idx_t]) and not np.isinf(Q_grid[idx_p,idx_t]) and not np.isinf(Bs_grid[idx_p,idx_t]) and not np.isinf(Bs_grid[idx_p,idx_t]) and not np.isinf(Be_grid[idx_p,idx_t]) and not np.isinf(Be_grid[idx_p,idx_t]) and 0.0>Q_grid[idx_p,idx_t]:
			Expansion_grid[idx_p,idx_t]=max(Bs_grid[idx_p,idx_t],Be_grid[idx_p,idx_t])/min(Bs_grid[idx_p,idx_t],Be_grid[idx_p,idx_t])


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig1=plt.figure(figsize=(10,5))
color_plot1=plt.pcolormesh(grid_phi*RAD2DEG,grid_theta*RAD2DEG-90.0,(np.log(Expansion_grid+0.1)/np.log(10.0))[:-1,:-1],rasterized=True,vmin=1.4,vmax=3.0,cmap="hot_r")
if plot_contour:
	contour1=plt.contour(grid_phi*RAD2DEG,grid_theta*RAD2DEG-90.0,np.log(Expansion_grid+0.1)/np.log(10.0),[1.0],colors=["blue"])

cbar1=fig1.colorbar(color_plot1)
cbar1.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar1.set_label(label=r"$\mathrm{log}_{10}(B_1/B_2)$",fontsize=20)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)

plt.ylim(theta_limits)
plt.xlim(phi_limits)


plt.savefig("Expansion_Factor.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)
plt.show()
