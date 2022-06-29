
"""
Calculate amount of reconnected magnetic flux 
Requires output fieldlines program (not included in this repository) and connectivity maps to be premade
"""

bfield_directory="./"
bfield_steps=[0000001,0000002,0000003,0000004,0000005]

connection_filename_base="./PFLS/Connection_0000000_"
fieldlines_filename_base="./PFLS/"


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from scipy.interpolate import griddata


connection_colors=[(171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)]
connection_cmap=matplotlib.colors.ListedColormap(((171/255,210/255,229/255),(244/255,166/255,131/255),(102/255,102/255,102/255),(143/255,145/255,0/255),(1.0,128/255,0)))


reconnected_flux_tot=np.zeros((5,len(bfield_steps)))
reconnected_flux_abs=np.zeros((5,len(bfield_steps)))
reconnected_flux_num=np.zeros((5,len(bfield_steps)),dtype=np.int32)
for idx_f in range(len(bfield_steps)):
	theta_grid,phi_grid,connection_map=load_connection_map(connection_filename_base+str(bfield_steps[idx_f])+".bin")
	#time,ntblks,nlblks,coord_logR,coord_theta,coord_phi,B_field=read_bfield_file(bfield_directory,"bfield."+str(bfield_steps[idx_f]))
	R,theta,phi,B_ends,flines,mapping,Q=read_fieldline_file(fieldlines_filename_base+str(bfield_steps[idx_f])+"/fieldline.r1l."+str(bfield_steps[idx_f]))
	delta_theta=theta_grid[0,1]-theta_grid[0,0]
	delta_phi=phi_grid[1,0]-phi_grid[0,0]
	for idx_p in range(len(phi_grid[:,0])):
		for idx_t in range(len(theta_grid[0,:])):
			if abs(mapping[idx_t,idx_p,0,0]-1.0)<0.1 and abs(theta[idx_t]-mapping[idx_t,idx_p,0,1])<1.0 and abs(phi[idx_p]-mapping[idx_t,idx_p,0,2])<1.0:
				idx_m=0
				idx_r=int(connection_map[idx_p,idx_t]//2)
				reconnected_flux_tot[idx_r,idx_f]+=B_ends[idx_t,idx_p,idx_m,0]*delta_theta*delta_phi*np.sin(theta_grid[0,idx_t])
				reconnected_flux_abs[idx_r,idx_f]+=abs(B_ends[idx_t,idx_p,idx_m,0])*delta_theta*delta_phi*np.sin(theta_grid[0,idx_t])
				reconnected_flux_num[idx_r,idx_f]+=1
			elif abs(mapping[idx_t,idx_p,1,0]-1.0)<0.1 and abs(theta[idx_t]-mapping[idx_t,idx_p,1,1])<1.0 and abs(phi[idx_p]-mapping[idx_t,idx_p,1,2])<1.0:
				idx_m=1
				idx_r=int(connection_map[idx_p,idx_t]//2)
				reconnected_flux_tot[idx_r,idx_f]+=B_ends[idx_t,idx_p,idx_m,0]*delta_theta*delta_phi*np.sin(theta_grid[0,idx_t])
				reconnected_flux_abs[idx_r,idx_f]+=abs(B_ends[idx_t,idx_p,idx_m,0])*delta_theta*delta_phi*np.sin(theta_grid[0,idx_t])
				reconnected_flux_num[idx_r,idx_f]+=1
			else:
				print("Error",idx_t,idx_p)
	print(reconnected_flux_tot[0,idx_f],reconnected_flux_tot[1,idx_f],reconnected_flux_tot[2,idx_f],reconnected_flux_tot[3,idx_f],reconnected_flux_tot[4,idx_f],reconnected_flux_abs[0,idx_f],reconnected_flux_abs[1,idx_f],reconnected_flux_abs[2,idx_f],reconnected_flux_abs[3,idx_f],reconnected_flux_abs[4,idx_f],reconnected_flux_num[0,idx_f],reconnected_flux_num[1,idx_f],reconnected_flux_num[2,idx_f],reconnected_flux_num[3,idx_f],reconnected_flux_num[4,idx_f])

plt.figure()
plt.plot(reconnected_flux_tot[2,:],color=connection_colors[2])
plt.plot(reconnected_flux_tot[3,:],color=connection_colors[3])

plt.figure()
plt.plot(reconnected_flux_abs[2,:],color=connection_colors[2])
plt.plot(reconnected_flux_abs[3,:],color=connection_colors[3])

plt.show()
