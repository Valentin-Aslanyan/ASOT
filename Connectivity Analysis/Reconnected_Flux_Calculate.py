
#TODO - wrong order on B_R_func!!
solar_Radius=1.0*7E10
border_px=10
filter_px=5
filter_thresh=5
R_start    =[+1.000,+1.000,+1.000,+1.000,+1.000,+1.000]
theta_start=[+1.998,+2.009,+1.994,+1.965,+1.971,+1.979]
phi_start  =[-0.129,-0.144,-0.163,-0.124,-0.138,-0.100]


frames_per_step=3
frames_per_sec=2
pad_start_frames=6
pad_end_frames=6


QSL_directories=[
"./PFLS/0046491/",
"./PFLS/0047589/",
"./PFLS/0048687/",
"./PFLS/0049786/",
"./PFLS/0050887/",
"./PFLS/0051993/"]



import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from mayavi import mlab
from subprocess import call
from scipy.interpolate import interp2d


num_timesteps=len(QSL_directories)
if filter_px % 2==0:
	filter_px+=1
flux_integral=np.zeros((num_timesteps-1))


R1,theta_grid1,phi_grid1,Q1=parse_QSL_Rbinfile(QSL_directories[0]+'qslR1.bin')
Q1sign_list=[np.sign(Q1)]
theta_grid1_filtered=theta_grid1[border_px:-border_px,border_px:-border_px]
phi_grid1_filtered=phi_grid1[border_px:-border_px,border_px:-border_px]
for idx in range(1,num_timesteps):
	R1,theta_grid1,phi_grid1,Q1=parse_QSL_Rbinfile(QSL_directories[idx]+'qslR1.bin')
	Q1sign_list.append(np.sign(Q1))
	reconnected=0.5-0.5*Q1sign_list[idx-1]*Q1sign_list[idx]
	reconnected_filtered=np.copy(reconnected)
	for idx_t in range(filter_px//2,len(reconnected[0,:])-filter_px):
		for idx_p in range(filter_px//2,len(reconnected[:,0])-filter_px):
			if sum(reconnected[idx_p-filter_px//2:idx_p+filter_px//2+1,idx_t-filter_px//2:idx_t+filter_px//2+1].flatten())<filter_thresh:
				reconnected_filtered[idx_p,idx_t]=0.0
	reconnected_filtered=reconnected_filtered[border_px:-border_px,border_px:-border_px]
	phi,theta,R,B=parse_QSL_folder(QSL_directories[idx])
	phi=(phi-180.0)*np.pi/180.0
	theta=(theta+90.0)*np.pi/180.0
	B_R_func=interp2d(phi,theta,B[0,:,:,0])

	for idx_t in range(len(reconnected_filtered[0,:])):
		for idx_p in range(len(reconnected_filtered[:,0])):
			if reconnected_filtered[idx_p,idx_t]==1.0:
				flux_integral[idx-1]+=solar_Radius*solar_Radius*np.sin(theta_grid1_filtered[0,idx_t])*(theta_grid1_filtered[0,idx_t]-theta_grid1_filtered[0,idx_t-1])*(phi_grid1_filtered[idx_p,0]-phi_grid1_filtered[idx_p-1,0])*B_R_func(phi_grid1_filtered[idx_p,0],theta_grid1_filtered[0,idx_t])
			

X_grid1=np.sin(theta_grid1)*np.cos(phi_grid1)
Y_grid1=np.sin(theta_grid1)*np.sin(phi_grid1)
Z_grid1=np.cos(theta_grid1)


plt.figure()
plt.plot(flux_integral)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\Phi$ [G cm$^2$]",fontsize=20)
plt.xlabel(r"$t$ [a.u.]",fontsize=20)


for idx in range(len(flux_integral)):
	print(str(flux_integral[idx])+", ",end="")
print()




plt.show()



