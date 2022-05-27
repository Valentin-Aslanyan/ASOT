

"""
Evaluate connectivity (e.g. interchange reconnection) based on surface flows
At photosphere, or traced up to higher radius
"""


use_Fortran_routines=True

solar_Radius=1.0*7E10
CNT_file="./arms.cnt"
start_file="./PFLS/0000000/qslR1.qsl"
end_file="./PFLS/1000000/qslR1.qsl"
connection_filename="./PFLS/Connection_0000000_1000000_r1.bin"
start_time=0.0
end_time=1000.0
delta_t=10.0

#Set the following if computing connectivity map away from photosphere
#To have a connection map at R!=1, must first generate one at R==1 and state it below
connection_radius=1.0	#At photosphere if =None or =1.0 or =False
custom_theta=[-20,50,2048]	#If above is true, put None or list: [lower,upper,number] (use degrees)
custom_phi=[-70,70,4096]	#e.g. [-30.0,30.0,2048]
bfield_step=2531467
generate_fieldlines_file=False
fieldlines_file="./PFLS/fieldline.1000000"
surface_connection_filename="./PFLS/Connection_0000000_1000000_r1.bin"


import sys
sys.path[:0]=['/Change/This/Path']
from ARMS_ASOT_Functions import *
from subprocess import call

if use_Fortran_routines:
	if connection_radius==None or connection_radius==False or connection_radius==1.0:
		Connectivity_CNT_file=open("Connectivity_Map_Generate.cnt","w")
		print("!arms.cnt file",file=Connectivity_CNT_file)
		print(CNT_file,file=Connectivity_CNT_file)
		print("!start Q file",file=Connectivity_CNT_file)
		print(start_file,file=Connectivity_CNT_file)
		print("!end Q file",file=Connectivity_CNT_file)
		print(end_file,file=Connectivity_CNT_file)
		print("!Output file",file=Connectivity_CNT_file)
		print(connection_filename,file=Connectivity_CNT_file)
		print("!start time",file=Connectivity_CNT_file)
		print(str(start_time),file=Connectivity_CNT_file)
		print("!end time",file=Connectivity_CNT_file)
		print(str(end_time),file=Connectivity_CNT_file)
		print("!delta_t",file=Connectivity_CNT_file)
		print(str(delta_t),file=Connectivity_CNT_file)
		print("!solar_Radius",file=Connectivity_CNT_file)
		print(str(solar_Radius),file=Connectivity_CNT_file)
		Connectivity_CNT_file.close()
		call_result=call(["./Connectivity_Map_Generate"])
	else:
		theta_grid,phi_grid,connection_map=load_connection_map(surface_connection_filename)
		R_lower=1.0
		theta_grid_l=min(theta_grid[0,:])
		theta_grid_r=max(theta_grid[0,:])
		theta_grid_delta=(theta_grid_r-theta_grid_l)/(len(theta_grid[0,:])-1)
		phi_grid_l=min(phi_grid[:,0])
		phi_grid_r=max(phi_grid[:,0])
		phi_grid_delta=(phi_grid_r-phi_grid_l)/(len(phi_grid[:,0])-1)
		if custom_theta==None or custom_phi==None:
			new_theta_grid=theta_grid[0,:]
			new_phi_grid=phi_grid[:,0]
		else:
			new_theta_grid=np.linspace((custom_theta[0]+90.0)*DEG2RAD,(custom_theta[1]+90.0)*DEG2RAD,num=custom_theta[2])
			new_phi_grid=np.linspace(custom_phi[0]*DEG2RAD,custom_phi[1]*DEG2RAD,num=custom_phi[2])

		if generate_fieldlines_file:
			fieldline_CNT_file=open("fieldline.cnt","w")
			print(" .true.		  	  !spherical?",file=fieldline_CNT_file)
			print(" .false.,.false.,.true.   !x,y,z periodic?",file=fieldline_CNT_file)
			print(" .false.			  !save field line positions?",file=fieldline_CNT_file)
			print(" 0, "+str(int(len(new_theta_grid)-1))+","+str(int(len(new_phi_grid)-1))+"           !n1,n2,n3",file=fieldline_CNT_file)
			print(" "+str(connection_radius)+", "+str(new_theta_grid[0]*RAD2DEG-90.0)+", "+str(new_phi_grid[0]*RAD2DEG)+"     !rlo",file=fieldline_CNT_file)
			print(" "+str(connection_radius)+", "+str(new_theta_grid[-1]*RAD2DEG-90.0)+", "+str(new_phi_grid[-1]*RAD2DEG)+"     !rhi",file=fieldline_CNT_file)
			print(" 1500, -0.02, 0.05, 1e-6    !nsteps,lstart,lmax, error bound",file=fieldline_CNT_file)
			print(" .false.                  !turnmap?",file=fieldline_CNT_file)
			print(" 1.0, 30.0	 	  !top,bot - ignored for spherical",file=fieldline_CNT_file)
			print(" 1                        !no. snapshots ",file=fieldline_CNT_file)
			print(str(int(bfield_step)),file=fieldline_CNT_file)
			fieldline_CNT_file.close()

			call_result=call(["./fieldline_asot"])
			call_result=call(["cp","./PFLS/fieldline."+f'{bfield_step:07}',fieldlines_file])

		R_fl,theta_fl,phi_fl,B_ends,flines,mapping,Q=read_fieldline_file(fieldlines_file)
		connection_map_new=np.zeros((len(phi_fl),len(theta_fl)),dtype=np.int32)
		for idx_p in range(len(phi_fl)):
			for idx_t in range(len(theta_fl)):
				if abs(1.0-mapping[idx_t,idx_p,0,0]/R_lower)<0.1:		#At least one end tied down	
					if abs(1.0-mapping[idx_t,idx_p,1,0]/R_lower)<0.1:	#Second end also tied down
						if (mapping[idx_t,idx_p,0,1]+90.0)*DEG2RAD>=theta_grid_l and (mapping[idx_t,idx_p,0,1]+90.0)*DEG2RAD<=theta_grid_r and mapping[idx_t,idx_p,0,2]*DEG2RAD>=phi_grid_l and mapping[idx_t,idx_p,0,2]*DEG2RAD<=phi_grid_r:
							idx_t2=int(np.round(((mapping[idx_t,idx_p,0,1]+90.0)*DEG2RAD-theta_grid_l)/theta_grid_delta))
							idx_p2=int(np.round((mapping[idx_t,idx_p,0,2]*DEG2RAD-phi_grid_l)/phi_grid_delta))
							type_temp1=connection_map[idx_p2,idx_t2]
						else:
							type_temp1=3
						if (mapping[idx_t,idx_p,1,1]+90.0)*DEG2RAD>=theta_grid_l and (mapping[idx_t,idx_p,1,1]+90.0)*DEG2RAD<=theta_grid_r and mapping[idx_t,idx_p,1,2]*DEG2RAD>=phi_grid_l and mapping[idx_t,idx_p,1,2]*DEG2RAD<=phi_grid_r:
							idx_t2=int(np.round(((mapping[idx_t,idx_p,1,1]+90.0)*DEG2RAD-theta_grid_l)/theta_grid_delta))
							idx_p2=int(np.round((mapping[idx_t,idx_p,1,2]*DEG2RAD-phi_grid_l)/phi_grid_delta))
							type_temp2=connection_map[idx_p2,idx_t2]
						else:
							type_temp2=3
						if type_temp1!=3:
							connection_map_new[idx_p,idx_t]=type_temp1
						elif type_temp2!=3:
							connection_map_new[idx_p,idx_t]=type_temp2
						else:
							connection_map_new[idx_p,idx_t]=3
					else: 							#Second end not tied down
						if (mapping[idx_t,idx_p,0,1]+90.0)*DEG2RAD>=theta_grid_l and (mapping[idx_t,idx_p,0,1]+90.0)*DEG2RAD<=theta_grid_r and mapping[idx_t,idx_p,0,2]*DEG2RAD>=phi_grid_l and mapping[idx_t,idx_p,0,2]*DEG2RAD<=phi_grid_r:
							idx_t2=int(np.round(((mapping[idx_t,idx_p,0,1]+90.0)*DEG2RAD-theta_grid_l)/theta_grid_delta))
							idx_p2=int(np.round((mapping[idx_t,idx_p,0,2]*DEG2RAD-phi_grid_l)/phi_grid_delta))
							connection_map_new[idx_p,idx_t]=connection_map[idx_p2,idx_t2]
						else:
							connection_map_new[idx_p,idx_t]=0
				elif abs(1.0-mapping[idx_t,idx_p,1,0]/R_lower)<0.1:				#Only one end tied down
					if (mapping[idx_t,idx_p,1,1]+90.0)*DEG2RAD>=theta_grid_l and (mapping[idx_t,idx_p,1,1]+90.0)*DEG2RAD<=theta_grid_r and mapping[idx_t,idx_p,1,2]*DEG2RAD>=phi_grid_l and mapping[idx_t,idx_p,1,2]*DEG2RAD<=phi_grid_r:
						idx_t2=int(np.round(((mapping[idx_t,idx_p,1,1]+90.0)*DEG2RAD-theta_grid_l)/theta_grid_delta))
						idx_p2=int(np.round((mapping[idx_t,idx_p,1,2]*DEG2RAD-phi_grid_l)/phi_grid_delta))
						connection_map_new[idx_p,idx_t]=connection_map[idx_p2,idx_t2]
					else:
						connection_map_new[idx_p,idx_t]=0
				else:						#Neither end tied down
					connection_map_new[idx_p,idx_t]=8

		outfile=open(connection_filename,'wb')
		outfile.write(struct.pack('i',len(phi_fl)))
		outfile.write(struct.pack('i',len(theta_fl)))
		for idx_p in range(len(phi_fl)):
			outfile.write(struct.pack('f',new_phi_grid[idx_p]))
		for idx_t in range(len(theta_fl)):
			outfile.write(struct.pack('f',new_theta_grid[idx_t]))
		for idx_p in range(len(phi_fl)):
			outfile.write(bytes(connection_map_new[idx_p,:]))
		outfile.close()
		
else:
	save_connection_map(connection_filename,CNT_file,start_file,end_file,start_time,end_time,delta_t,solar_Radius)




