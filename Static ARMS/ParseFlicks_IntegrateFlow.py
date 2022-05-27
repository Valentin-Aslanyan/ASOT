

target_R=1.0	#degrees
file_directory="./"
flicks_files=['flicks.0044295',
 'flicks.0045393',
 'flicks.0046491',
 'flicks.0047589',
 'flicks.0048687',
 'flicks.0049786',
 'flicks.0050889',
 'flicks.0052001',
 'flicks.0053135',
 'flicks.0054317',
 'flicks.0055578',
 'flicks.0056909',
 'flicks.0058389',
 'flicks.0060201',
 'flicks.0062329',
 'flicks.0064701',
 'flicks.0067271',
 'flicks.0070167',
 'flicks.0073541',
 'flicks.0077304',
 'flicks.0081418',
 'flicks.0086001',
 'flicks.0091110',
 'flicks.0096469',
 'flicks.0102673',
 'flicks.0108951',
 'flicks.0115218',
 'flicks.0121501',
 'flicks.0127787',
 'flicks.0134087',
 'flicks.0140389',
 'flicks.0146688',
 'flicks.0152990',
 'flicks.0159295',
 'flicks.0165602',
 'flicks.0171913',
 'flicks.0178228',
 'flicks.0184543',
 'flicks.0190857',
 'flicks.0197169',
 'flicks.0203482',
 'flicks.0209795',
 'flicks.0216109',
 'flicks.0222422',
 'flicks.0228736',
 'flicks.0235046',
 'flicks.0241352',
 'flicks.0247660',
 'flicks.0253968',
 'flicks.0260279',
 'flicks.0266591',
 'flicks.0272901',
 'flicks.0279208',
 'flicks.0285517',
 'flicks.0291824',
 'flicks.0298128',
 'flicks.0304420',
 'flicks.0310709',
 'flicks.0317002',
 'flicks.0323299',
 'flicks.0329595',
 'flicks.0335888',
 'flicks.0342185',
 'flicks.0348482',
 'flicks.0354780',
 'flicks.0361077',
 'flicks.0367380',
 'flicks.0373687',
 'flicks.0379993',
 'flicks.0386297',
 'flicks.0392598',
 'flicks.0398897',
 'flicks.0405194',
 'flicks.0411489',
 'flicks.0417783',
 'flicks.0424077',
 'flicks.0430366',
 'flicks.0436654',
 'flicks.0442948',
 'flicks.0449245',
 'flicks.0455538',
 'flicks.0461832',
 'flicks.0468132',
 'flicks.0474434',
 'flicks.0480737',
 'flicks.0487035',
 'flicks.0493325',
 'flicks.0499615',
 'flicks.0505909',
 'flicks.0512207',
 'flicks.0518501',
 'flicks.0524796',
 'flicks.0531090',
 'flicks.0537385',
 'flicks.0543684',
 'flicks.0549988']


import sys
sys.path[:0]=['/Change/This/Path']
from ASOT_Functions_Python import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,flicks_files[0])
new_theta,new_phi,new_data=R_slice(target_R,coord_logR,coord_theta,coord_phi,data)
v_integrated=new_data[:,1:4]
v2=new_data[:,1]**2+new_data[:,2]**2+new_data[:,3]**2
max_v=[max(v2.flatten())]

for fname in flicks_files[1:]:
	time,ntblks,nlblks,newgrd,coord_logR,coord_theta,coord_phi,data=read_flicks_file(file_directory,fname)
	new_theta,new_phi,new_data=R_slice(target_R,coord_logR,coord_theta,coord_phi,data)
	v_integrated+=new_data[:,1:4]
	v2=new_data[:,1]**2+new_data[:,2]**2+new_data[:,3]**2
	max_v.append(max(v2.flatten()))


new_v=np.sqrt(v_integrated[:,0]**2+v_integrated[:,1]**2+v_integrated[:,2]**2)

min_new_theta=min(new_theta)
max_new_theta=max(new_theta)
min_new_phi=min(new_phi)
max_new_phi=max(new_phi)

from scipy.interpolate import griddata
grid_theta=np.linspace(min_new_theta,max_new_theta,num=500)
grid_phi=np.linspace(min_new_phi,max_new_phi,num=400)
grid_theta,grid_phi=np.meshgrid(grid_theta,grid_phi)
grid_theta2=np.linspace(110/180.0*np.pi,145/180.0*np.pi,num=500)
grid_phi2=np.linspace(-40/180.0*np.pi,60/180.0*np.pi,num=400)
grid_theta2,grid_phi2=np.meshgrid(grid_theta2,grid_phi2)
v_grid=griddata((new_theta,new_phi), new_v, (grid_theta2,grid_phi2), method='linear')

fig2=plt.figure("pcolormesh v",figsize=(9,4))
ax2=fig2.gca()
plt.title(r"$|\sum v|,\quad$ $R="+"{:.1f}".format(target_R)+"$",fontsize=20)
colour_plot2=plt.pcolormesh(grid_phi2/np.pi*180.0,grid_theta2/np.pi*180.0,v_grid,cmap='Reds',rasterized=True)
plt.tick_params(axis='both', which='major',labelsize=19,direction='in',bottom=True, top=True, left=True, right=True)
plt.ylabel(r"$\theta$ [$^{\circ}$]",fontsize=20)
plt.xlabel(r"$\phi$ [$^{\circ}$]",fontsize=20)
cbar2=fig2.colorbar(colour_plot2)#,ticks=[-4,-3,-2,-1,0])
cbar2.ax.tick_params(labelsize=19,direction='in', left=True, right=True)
cbar2.set_label(label=r"$|\sum v|$",fontsize=20)
#plt.axis('off')
#plt.xlim([0,3.1])
#plt.ylim([-3,3])
#plt.savefig("IntegralV.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

fig1=plt.figure("Max |v|",figsize=(9,4))
ax1=fig1.gca()
plt.plot(np.sqrt(np.array(max_v)))

plt.show()



