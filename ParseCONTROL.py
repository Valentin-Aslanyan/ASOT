

target_profile_idx=[0,1,2]	#must be list
base_directory="./"


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

if len(sys.argv)<2:
	print("Specify subfolder as command line argument")
else:
	target_directory=os.path.join(base_directory,sys.argv[1])
	CNT_file=os.path.join(target_directory,'arms.cnt')

	Time_profiles,parameters=parse_sflow_CNT(CNT_file)


	fig=plt.figure("1",figsize=(10,6))
	plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
	plt.ylabel(r"Time amplitude",fontsize=20)
	plt.xlabel(r"$t$ [s]",fontsize=20)


	time_integral=0.0
	for idx in target_profile_idx:
		if np.shape(parameters)[0]>idx:
			t_p=2.0*np.pi*parameters[idx,18]/(parameters[idx,16]-parameters[idx,15])
			t=np.linspace(parameters[idx,15],parameters[idx,16],num=2000)
			if Time_profiles[idx]==0:
				t_func=np.ones((len(t)))
			if Time_profiles[idx]==1:
				t_func=(t-parameters[idx,17])/(parameters[idx,16]-parameters[idx,15])
			if Time_profiles[idx]==2:
				t_func=0.5*(1.0-np.cos(t_p*(t-parameters[idx,17])))
			if Time_profiles[idx]==3:
				t_func=np.sin(t_p*(t-parameters[idx,17]))
			plt.plot(t,t_func)
			time_integral+=np.trapz(t_func,x=t)
	print("Time integral [trapezoidal]: ",time_integral)

	plt.savefig("Time_Trace.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)
	plt.show()


