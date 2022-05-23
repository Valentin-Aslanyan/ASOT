
"""
Attempt to predict grid shape from variables in arms.cnt without running ARMS itself
"""


plot_defined_region=False #Red line showing defining regions
target_R=1.0
target_phi=1.0


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#                 rleft     rrigh     tleft     trigh     pleft     prigh    grdscl
grid_params=[     1.0,       3.0,    0.0625,   0.9375,     -1.0,    +1.0,    7.0e+10]

#            nblockr  nblockt  nblockp
init_data=[       4,        9,      20]

#                 dmv1l    dmv1r    dmv2l   dmv2r    dmv3l    dmv3r  lrefdmv
refine_data=[[7.0e+10,  2.1e+11,    0.2,    +0.7,     -0.45,      0.5,        2],
   [7.0e+10,  9.1e+10,   0.15,    +0.35,    -0.25,      0.1,        4],
   [7.0e+10,  9.1e+10,   0.15,    +0.4,       0.0,      0.4,        4],
   [7.0e+10, 7.026e+10, 0.4,    +0.5,      -0.15,      0.1,        5],
   [7.0e+10,  7.8e+10, 0.13,    +0.42,     -0.3,      0.4,       5],
   [7.8e+10,  1.0e+11, 0.15,    +0.4,     -0.3,      0.4,       5]]


#Preliminaries
max_ref_level=1
refine_data_proc=[]
for idx_ref in range(len(refine_data)):
	temp=refine_data[idx_ref].copy()
	temp[0]=np.log(temp[0]/grid_params[6])-np.log(grid_params[0])
	temp[1]=np.log(temp[1]/grid_params[6])-np.log(grid_params[0])
	temp[2]=temp[2]-grid_params[2]
	temp[3]=temp[3]-grid_params[2]
	temp[4]=temp[4]-grid_params[4]
	temp[5]=temp[5]-grid_params[4]
	refine_data_proc.append(temp)
	if max_ref_level<refine_data[idx_ref][6]:
		max_ref_level=refine_data[idx_ref][6]



#Basic grid, use units of Minimum Spacing (i.e. 1MS = size of most refined block)
r_width=np.log(grid_params[1])-np.log(grid_params[0])
t_width=grid_params[3]-grid_params[2]
p_width=grid_params[5]-grid_params[4]
r_spacing=r_width/init_data[0]/2**(max_ref_level-1)
t_spacing=t_width/init_data[1]/2**(max_ref_level-1)
p_spacing=p_width/init_data[2]/2**(max_ref_level-1)
max_r=2**(max_ref_level-1)*init_data[0]
max_t=2**(max_ref_level-1)*init_data[1]
max_p=2**(max_ref_level-1)*init_data[2]
coord_r=[]
coord_t=[]
coord_p=[]
ref_level=[]
for idx_r in range(init_data[0]):
	for idx_t in range(init_data[1]):
		for idx_p in range(init_data[2]):
			step_size=2**(max_ref_level-1)
			coord_r.append([idx_r*step_size,(idx_r+1)*step_size])
			coord_t.append([idx_t*step_size,(idx_t+1)*step_size])
			coord_p.append([idx_p*step_size,(idx_p+1)*step_size])
			ref_level.append(1)

#Adjust refinement specs to MS units; pad spec to next (up or down) block of that size
for idx_ref in range(len(refine_data)):
	current_spacing=2**(max_ref_level+1-refine_data_proc[idx_ref][6])
	refine_data_proc[idx_ref][0]=int(max(0,np.floor(refine_data_proc[idx_ref][0]/r_spacing/current_spacing)*current_spacing))
	refine_data_proc[idx_ref][1]=int(min(max_r,np.ceil(refine_data_proc[idx_ref][1]/r_spacing/current_spacing)*current_spacing))
	refine_data_proc[idx_ref][2]=int(max(0,np.floor(refine_data_proc[idx_ref][2]/t_spacing/current_spacing)*current_spacing))
	refine_data_proc[idx_ref][3]=int(min(max_t,np.ceil(refine_data_proc[idx_ref][3]/t_spacing/current_spacing)*current_spacing))
	refine_data_proc[idx_ref][4]=int(max(0,np.floor(refine_data_proc[idx_ref][4]/p_spacing/current_spacing)*current_spacing))
	refine_data_proc[idx_ref][5]=int(min(max_p,np.ceil(refine_data_proc[idx_ref][5]/p_spacing/current_spacing)*current_spacing))

#Refinement
for idx_l in range(2,max_ref_level+1):
	new_spacing=2**(max_ref_level-idx_l)
	for idx_ref in range(len(refine_data_proc)):
		new_ref_level=refine_data_proc[idx_ref][6]
		if new_ref_level>=idx_l:
			#Adjust r to buffer adjacent blocks
			r_ref_low=refine_data_proc[idx_ref][0]
			r_ref_high=refine_data_proc[idx_ref][1]
			for idx3 in range(new_ref_level-1,idx_l-1,-1):
				current_spacing=2**(max_ref_level-idx3)
				if r_ref_low % (current_spacing*2)==0:
					r_ref_low-=2*current_spacing
				else:
					r_ref_low-=1*current_spacing
				if r_ref_high % (current_spacing*2)==0:
					r_ref_high+=2*current_spacing
				else:
					r_ref_high+=1*current_spacing
			r_ref_low=max(0,r_ref_low)
			r_ref_high=min(max_r,r_ref_high)

			#Adjust t to buffer adjacent blocks
			t_ref_low=refine_data_proc[idx_ref][2]
			t_ref_high=refine_data_proc[idx_ref][3]
			for idx3 in range(new_ref_level-1,idx_l-1,-1):
				current_spacing=2**(max_ref_level-idx3)
				if t_ref_low % (current_spacing*2)==0:
					t_ref_low-=2*current_spacing
				else:
					t_ref_low-=1*current_spacing
				if t_ref_high % (current_spacing*2)==0:
					t_ref_high+=2*current_spacing
				else:
					t_ref_high+=1*current_spacing
			t_ref_low=max(0,t_ref_low)
			t_ref_high=min(max_t,t_ref_high)

			#Adjust p to buffer adjacent blocks
			p_ref_low=refine_data_proc[idx_ref][4]
			p_ref_high=refine_data_proc[idx_ref][5]
			for idx3 in range(new_ref_level-1,idx_l-1,-1):
				current_spacing=2**(max_ref_level-idx3)
				if p_ref_low % (current_spacing*2)==0:
					p_ref_low-=2*current_spacing
				else:
					p_ref_low-=1*current_spacing
				if p_ref_high % (current_spacing*2)==0:
					p_ref_high+=2*current_spacing
				else:
					p_ref_high+=1*current_spacing
			p_ref_low=max(0,p_ref_low)
			p_ref_high=min(max_p,p_ref_high)

			for idx in range(len(coord_r)-1,-1,-1):
				if ref_level[idx]<idx_l and r_ref_low<=coord_r[idx][0] and r_ref_high>=coord_r[idx][1] and t_ref_low<=coord_t[idx][0] and t_ref_high>=coord_t[idx][1] and p_ref_low<=coord_p[idx][0] and p_ref_high>=coord_p[idx][1]: 
					r_start=coord_r[idx][0]
					t_start=coord_t[idx][0]
					p_start=coord_p[idx][0]
					coord_r.pop(idx)
					coord_t.pop(idx)
					coord_p.pop(idx)
					ref_level.pop(idx)
					for idx_r in range(2):
						for idx_t in range(2):
							for idx_p in range(2):
								r_low=r_start+idx_r*new_spacing
								r_high=r_start+(idx_r+1)*new_spacing
								t_low=t_start+idx_t*new_spacing
								t_high=t_start+(idx_t+1)*new_spacing
								p_low=p_start+idx_p*new_spacing
								p_high=p_start+(idx_p+1)*new_spacing
								coord_r.append([r_low,r_high])
								coord_t.append([t_low,t_high])
								coord_p.append([p_low,p_high])
								ref_level.append(idx_l)



print("Num blocks: ",len(coord_r))


#Convert to real units
coord_logR=[]
coord_theta=[]
coord_phi=[]
for idx in range(len(coord_r)):
	coord_logR.append([coord_r[idx][0]*r_spacing+np.log(grid_params[0]),coord_r[idx][1]*r_spacing+np.log(grid_params[0])])
	coord_theta.append([(1.0-(coord_t[idx][1]*t_spacing+grid_params[2]))*np.pi,(1.0-(coord_t[idx][0]*t_spacing+grid_params[2]))*np.pi])
	coord_phi.append([(coord_p[idx][0]*p_spacing+grid_params[4])*np.pi,(coord_p[idx][1]*p_spacing+grid_params[4])*np.pi])


plt.figure("Phi slice",figsize=(3.1*1.5,6*1.5))
for idx in range(len(coord_logR)):
	if coord_phi[idx][0]<=target_phi/180.0*np.pi and coord_phi[idx][1]>=target_phi/180.0*np.pi:
		x1=np.exp(coord_logR[idx][0])*np.sin(coord_theta[idx][0])
		z1=np.exp(coord_logR[idx][0])*np.cos(np.pi-coord_theta[idx][0])
		x2=np.exp(coord_logR[idx][1])*np.sin(coord_theta[idx][0])
		z2=np.exp(coord_logR[idx][1])*np.cos(np.pi-coord_theta[idx][0])
		x3=np.exp(coord_logR[idx][1])*np.sin(coord_theta[idx][1])
		z3=np.exp(coord_logR[idx][1])*np.cos(np.pi-coord_theta[idx][1])
		x4=np.exp(coord_logR[idx][0])*np.sin(coord_theta[idx][1])
		z4=np.exp(coord_logR[idx][0])*np.cos(np.pi-coord_theta[idx][1])
		plt.plot([x1,x2],[z1,z2],color="grey",linewidth=1)
		plt.plot([x3,x4],[z3,z4],color="grey",linewidth=1)
		anglearc=np.linspace(coord_theta[idx][0],coord_theta[idx][1],num=100)
		arc1x=np.exp(coord_logR[idx][0])*np.sin(anglearc)
		arc1y=np.exp(coord_logR[idx][0])*np.cos(np.pi-anglearc)
		arc2x=np.exp(coord_logR[idx][1])*np.sin(anglearc)
		arc2y=np.exp(coord_logR[idx][1])*np.cos(np.pi-anglearc)
		plt.plot(arc1x,arc1y,color="grey",linewidth=1)
		plt.plot(arc2x,arc2y,color="grey",linewidth=1)
plt.axis('off')
plt.xlim([0,3.1])
plt.ylim([-3,3])
plt.savefig("Predicted_Grid_XZ.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)


plt.figure("R slice",figsize=(20,10))
for idx in range(len(coord_logR)):
	if np.exp(coord_logR[idx][0])<=target_R and np.exp(coord_logR[idx][1])>=target_R:
		th1=coord_theta[idx][0]/np.pi*180.0
		ph1=coord_phi[idx][0]/np.pi*180.0
		th2=coord_theta[idx][1]/np.pi*180.0
		ph2=coord_phi[idx][1]/np.pi*180.0
		plt.plot([ph1,ph1,ph2,ph2,ph1],[th1,th2,th2,th1,th1],color="grey",linewidth=1)

if plot_defined_region:
	for idx_ref in range(len(refine_data)):
		if refine_data[idx_ref][0]/grid_params[6]<=target_R and refine_data[idx_ref][1]/grid_params[6]>=target_R:
			th1=(1.0-refine_data[idx_ref][2])*180.0
			ph1=refine_data[idx_ref][4]*180.0
			th2=(1.0-refine_data[idx_ref][3])*180.0
			ph2=refine_data[idx_ref][5]*180.0
			plt.plot([ph1,ph1,ph2,ph2,ph1],[th1,th2,th2,th1,th1],color="red",linewidth=1)

	"""
	for idx_ref in range(len(refine_data)):
		th1=(1.0-refine_data_proc[idx_ref][2]*t_spacing-grid_params[2])*180.0
		ph1=(refine_data_proc[idx_ref][4]*p_spacing+grid_params[4])*180.0
		th2=(1.0-refine_data_proc[idx_ref][3]*t_spacing-grid_params[2])*180.0
		ph2=(refine_data_proc[idx_ref][5]*p_spacing+grid_params[4])*180.0
		plt.plot([ph1,ph1,ph2,ph2,ph1],[th1,th2,th2,th1,th1],color="green",linewidth=1)

	for idx_ref in range(len(refine_data)):
		th1=(1.0-t_ref_low*t_spacing-grid_params[2])*180.0
		ph1=(p_ref_low*p_spacing+grid_params[4])*180.0
		th2=(1.0-t_ref_high*t_spacing-grid_params[2])*180.0
		ph2=(p_ref_high*p_spacing+grid_params[4])*180.0
		plt.plot([ph1,ph1,ph2,ph2,ph1],[th1,th2,th2,th1,th1],color="blue",linewidth=1)
	"""

plt.axis('off')
plt.savefig("Predicted_Grid_PhiTheta.pdf", format="pdf", dpi=100,bbox_inches='tight',pad_inches=0.1)

plt.show()


