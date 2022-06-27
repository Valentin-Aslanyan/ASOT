
"""
Display squashing factor from QSL Squasher (not included in this distribution)
"""

filename="./PFLS/0000000/qslR1.bin"
plot_PIL=False	#Polarity Inversion Line, at outer boundary


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


R_actual,theta_grid,phi_grid,Q=parse_QSL_Rbinfile(filename)

if plot_PIL:
	idx=len(filename)-1
	while idx>0 and filename[idx]!='/' and filename[idx]!="\\":
		idx-=1
	phi,theta,R,B=parse_QSL_folder(filename[:idx])
	phi=(phi-180.0)*DEG2RAD
	theta=(theta+90.0)*DEG2RAD
	theta_PIL,phi_PIL=get_Polarity_Inversion_Line(phi,theta,B[-2,:,:,0])

Q_grid=np.sign(Q)*np.log(abs(Q))/np.log(10.0)
Q_grid[np.isinf(Q_grid)]=np.nan
X_grid=np.sin(theta_grid)*np.cos(phi_grid)
Y_grid=np.sin(theta_grid)*np.sin(phi_grid)
Z_grid=np.cos(theta_grid)


fig=plt.figure("Q 2D",figsize=(20,10))
ax=fig.gca()
plt.title(r"$R=1R_{\odot}$",fontsize=20)
color_plot=plt.pcolormesh(phi_grid*RAD2DEG,theta_grid*RAD2DEG-90.0,Q_grid,cmap='RdBu_r',vmin=-5,vmax=5,rasterized=True)
if plot_PIL:
	plt.plot(phi_PIL*RAD2DEG,theta_PIL*RAD2DEG-90.0,color="black",linewidth=2)
	plt.plot(phi_PIL*RAD2DEG,theta_PIL*RAD2DEG-90.0,color="orange",linewidth=1)


plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar=fig.colorbar(color_plot)#,ticks=[-4,-3,-2,-1,0])
cbar.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar.set_label(label=r"$\mathrm{slog}(Q)$",fontsize=20)

plt.savefig("QSL.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.show()



