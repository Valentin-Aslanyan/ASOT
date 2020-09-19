

target_file="./PFLS/fieldline.r1.2112857"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


R,theta,phi,B_ends,flines,mapping,Q=read_fieldline_file(target_file)
phi_grid,theta_grid=np.meshgrid(phi,theta)
max_Q=max(abs(Q).flatten())
Q_grid=np.sign(Q)*np.log(np.clip(abs(Q),2.0,max_Q))/np.log(10.0)
Q_grid[np.isinf(Q_grid)]=np.nan

fig=plt.figure("Q 2D",figsize=(20,10))
ax=fig.gca()
#plt.title(r"$R=1R_{\odot}$",fontsize=20)

color_plot=plt.pcolormesh(phi_grid,theta_grid,Q_grid,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)

plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)

plt.savefig("Q_Fieldlines.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.show()




