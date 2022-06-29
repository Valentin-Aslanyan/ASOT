

filename1="./PFLS/0044295/qslR1.bin"
filename2="./PFLS/0233915/qslR1.bin"


theta_start=[+1.998,+2.009,+1.994,+1.965,+1.971,+1.979]
phi_start  =[-0.129,-0.144,-0.163,-0.124,-0.138,-0.100]

R_start    =+1.000
theta_start=+1.971
phi_start  =-0.138


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


#Field lines
phi1,theta1,R1,B1=parse_QSL_folder(filename1.replace("qslR1.bin",""))
phi1=(phi1-180.0)*np.pi/180.0
theta1=(theta1+90.0)*np.pi/180.0

phi2,theta2,R2,B2=parse_QSL_folder(filename2.replace("qslR1.bin",""))
phi2=(phi2-180.0)*np.pi/180.0
theta2=(theta2+90.0)*np.pi/180.0


field_line_start=np.array([R_start,theta_start,phi_start])
field_line_sph=field_line_spherical(field_line_start,R1,theta1,phi1,B1,1.0,2.9,step_size=1E-2)
field_line_X1=field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2])
field_line_Y1=-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2])
field_line_Z1=field_line_sph[:,0]*np.cos(field_line_sph[:,1])
field_line_sph=field_line_spherical(field_line_start,R2,theta2,phi2,B2,1.0,2.9,step_size=1E-2)
field_line_X2=field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.cos(field_line_sph[:,2])
field_line_Y2=-field_line_sph[:,0]*np.sin(field_line_sph[:,1])*np.sin(field_line_sph[:,2])
field_line_Z2=field_line_sph[:,0]*np.cos(field_line_sph[:,1])

delta_line=np.zeros((min(len(field_line_X1),len(field_line_X2))))
for idx in range(len(delta_line)):
	delta_line[idx]=np.sqrt((field_line_X1[idx]-field_line_X2[idx])**2+(field_line_Y1[idx]-field_line_Y2[idx])**2+(field_line_Z1[idx]-field_line_Z2[idx])**2)


#Q map
R1,theta1,phi1,Q1=parse_QSL_Rbinfile(filename1)
Q_grid1=np.sign(Q1)*np.log(abs(Q1))
Q_grid1[np.isinf(Q_grid1)]=np.nan
R2,theta2,phi2,Q2=parse_QSL_Rbinfile(filename2)
Q_grid2=np.sign(Q2)*np.log(abs(Q2))
Q_grid2[np.isinf(Q_grid2)]=np.nan

Q_signchange=np.zeros(np.shape(Q_grid1))
for idx_x in range(np.shape(Q_grid1)[0]):
	for idx_y in range(np.shape(Q_grid1)[1]):
		if not np.isnan(Q_grid1[idx_x,idx_y]) and not np.isnan(Q_grid2[idx_x,idx_y]) and Q_grid1[idx_x,idx_y]*Q_grid2[idx_x,idx_y]<0.0:
			Q_signchange[idx_x,idx_y]=1.0

Q_delta=Q_grid1[:,:]-Q_grid2[:,:]


plt.figure(figsize=(10,10))
plt.subplot(2,1,1)
plt.pcolormesh(phi1*180.0/np.pi,theta1*180.0/np.pi,Q_grid1,cmap='RdBu_r',vmin=-10,vmax=10)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.subplot(2,1,2)
plt.pcolormesh(phi2*180.0/np.pi,theta2*180.0/np.pi,Q_grid2,cmap='RdBu_r',vmin=-10,vmax=10)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.tight_layout()

plt.figure(figsize=(10,10))
plt.subplot(2,1,1)
plt.pcolormesh(phi1*180.0/np.pi,theta1*180.0/np.pi,Q_delta,cmap='PiYG_r',vmin=-10,vmax=10)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.subplot(2,1,2)
plt.pcolormesh(phi2*180.0/np.pi,theta2*180.0/np.pi,Q_signchange,cmap='binary',vmin=0,vmax=1)
plt.xlabel(r'$\phi$ [$^{\circ}$]',fontsize=20)
plt.ylabel(r'$\theta$ [$^{\circ}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.tight_layout()

plt.figure("r",figsize=(6,6))
plt.plot(np.sqrt(field_line_X1**2+field_line_Y1**2+field_line_Z1**2),color="red")
plt.plot(np.sqrt(field_line_X2**2+field_line_Y2**2+field_line_Z2**2),color="blue")
plt.xlabel(r'Path length [a.u.]',fontsize=20)
plt.ylabel(r'$R$ [$R_{\odot}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)


plt.figure("delta",figsize=(6,6))
plt.plot(delta_line,color="purple")
plt.xlabel(r'Path length [a.u.]',fontsize=20)
plt.ylabel(r'$|\vec{r}_1-\vec{r}_2|$ [$R_{\odot}$]',fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=18,direction="in",bottom=True, top=True, left=True, right=True)

plt.show()







