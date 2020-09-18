
filename1="./PFLS/0139195/qslR1.bin"
filename2="./PFLS/0139195/qslR3.bin"


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


R1,theta1,phi1,Q1=parse_QSL_Rbinfile(filename1)
Q_grid1=np.sign(Q1)*np.log(abs(Q1))
Q_grid1[np.isinf(Q_grid1)]=np.nan
R2,theta2,phi2,Q2=parse_QSL_Rbinfile(filename2)
Q_grid2=np.sign(Q2)*np.log(abs(Q2))
Q_grid2[np.isinf(Q_grid2)]=np.nan


plt.figure(figsize=(10,10))
plt.subplot(2,1,1)
plt.pcolormesh(phi1*180.0/np.pi,theta1*180.0/np.pi,Q_grid1,cmap='RdBu_r',vmin=-10,vmax=10)
plt.title(r"$R="+"{:.2f}".format(R1)+r"$ Mm",fontsize=20)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.subplot(2,1,2)
plt.pcolormesh(phi2*180.0/np.pi,theta2*180.0/np.pi,Q_grid2,cmap='RdBu_r',vmin=-10,vmax=10)
plt.title(r"$R="+"{:.2f}".format(R2)+r"$ Mm",fontsize=20)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.tight_layout()
plt.show()







