
"""
Plot squashing factor from QSL Squasher (not included in this distribution) at different times/radii alongside
"""

filename1="./PFLS/0000000/qslR1.bin"
filename2="./PFLS/0000000/qslR3.bin"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


R1,theta1,phi1,Q1=parse_QSL_Rbinfile(filename1)
Q_grid1=np.sign(Q1)*np.log(abs(Q1))/np.log(10.0)
Q_grid1[np.isinf(Q_grid1)]=np.nan
R2,theta2,phi2,Q2=parse_QSL_Rbinfile(filename2)
Q_grid2=np.sign(Q2)*np.log(abs(Q2))/np.log(10.0)
Q_grid2[np.isinf(Q_grid2)]=np.nan


plt.figure(figsize=(10,10))
plt.subplot(2,1,1)
plt.pcolormesh(phi1*RAD2DEG,theta1*RAD2DEG-90.0,Q_grid1,cmap='RdBu_r',vmin=-5,vmax=5)
plt.title(r"$R="+"{:.2f}".format(R1)+r"$ Mm",fontsize=20)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.subplot(2,1,2)
plt.pcolormesh(phi2*RAD2DEG,theta2*RAD2DEG-90.0,Q_grid2,cmap='RdBu_r',vmin=-5,vmax=5)
plt.title(r"$R="+"{:.2f}".format(R2)+r"$ Mm",fontsize=20)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.tight_layout()
plt.show()







